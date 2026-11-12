"""All 230 space groups: smallest test cell, constants, orientation and table coverage.

For every space group a primitive cell with Si and Ge on two general positions
is built (so the crystal has exactly that space group, not a higher one), in
spglib's standard setting. Then:

* auto finds the right number of independent constants for the Laue class,
  also after an arbitrary rotation, with well-conditioned deformations;
* the named function for that space group can represent the crystal's
  stiffness matrix, and the symmetry check warns exactly when it cannot.
"""
import warnings

import pytest
import spglib
from ase import Atoms
from ase.spacegroup import crystal
from scipy.spatial.transform import Rotation

from elastic_vasp import auto, orientation, poscar
from elastic_vasp.systems import SYSTEMS

GENERAL = [(0.1137, 0.2213, 0.3471), (0.3719, 0.0613, 0.1789)]

# (last space group of the range, number of constants, named system to use)
RANGES = [(2, 21, "triclinic"), (15, 13, "monoclinic"), (74, 9, "orthogonal"),
          (88, 7, "tetragonal_2"), (142, 6, "tetragonal_1"), (148, 7, "trigonal_2"),
          (167, 6, "trigonal_1"), (194, 5, "hexagonal"), (230, 3, "cubic")]
# trigonal I groups whose 2-fold axes are perpendicular to a in the standard
# setting ("31m" type): C15 instead of C14, so trigonal_2 is the named choice
TRIGONAL_31M = {149, 151, 153, 157, 159, 162, 163}


def expected(n):
    for last, count, named in RANGES:
        if n <= last:
            if named == "trigonal_1" and n in TRIGONAL_31M:
                named = "trigonal_2"
            return count, named
    raise ValueError(n)


def _cellpar(n):
    if n <= 2:
        return [4.1, 4.7, 5.3, 77, 83, 71]
    if n <= 15:
        return [4.1, 4.7, 5.3, 90, 101, 90]
    if n <= 74:
        return [4.1, 4.7, 5.3, 90, 90, 90]
    if n <= 142:
        return [4.1, 4.1, 5.3, 90, 90, 90]
    if n <= 194:
        return [4.1, 4.1, 5.3, 90, 90, 120]
    return [4.7, 4.7, 4.7, 90, 90, 90]


def smallest_cell(n):
    """Primitive cell of space group n (standard setting) with Si, Ge on general positions."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        conventional = crystal(["Si", "Ge"], basis=GENERAL, spacegroup=n,
                               cellpar=_cellpar(n), symprec=1e-4)
    cell = (conventional.cell.array, conventional.get_scaled_positions(), conventional.numbers)
    lattice, positions, numbers = spglib.standardize_cell(cell, to_primitive=True,
                                                          no_idealize=False, symprec=1e-4)
    return Atoms(numbers=numbers, cell=lattice, scaled_positions=positions, pbc=True)


@pytest.mark.parametrize("n", range(1, 231))
def test_space_group(n, tmp_path):
    atoms = smallest_cell(n)
    dataset = spglib.get_symmetry_dataset(
        (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers), symprec=1e-4)
    number = dataset.number if hasattr(dataset, "number") else dataset["number"]
    assert number == n  # exactly this space group, not a higher one
    count, named = expected(n)

    # auto: right number of constants, in the standard and in a random orientation
    spec, basis = auto.from_structure(atoms)
    assert len(spec.constants) == count
    rotated = atoms.copy()
    rotated.set_cell(rotated.cell.array @ Rotation.random(random_state=n).as_matrix().T,
                     scale_atoms=True)
    rotated_spec, _ = auto.from_structure(rotated)
    assert len(rotated_spec.constants) == count
    assert auto.amplification(rotated_spec) < 8

    # a rotated supercell: spglib alone sees only the operations that also map the
    # supercell onto itself, so this checks that the full point group is still used
    supercell = atoms.repeat((2, 1, 1))
    supercell.set_cell(supercell.cell.array @ Rotation.random(random_state=1000 + n).as_matrix().T,
                       scale_atoms=True)
    assert len(auto.from_structure(supercell)[0].constants) == count

    # the named function for this space group covers it, silently
    path = tmp_path / "POSCAR"
    poscar.write(path, atoms, "SG {}".format(n))
    assert orientation.table_covers(SYSTEMS[named], basis)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert orientation.check_table(named, str(path))

    # and the symmetry check warns for every named table that cannot represent it
    for name, table in SYSTEMS.items():
        if not orientation.table_covers(table, basis):
            with pytest.warns(UserWarning, match="cannot measure this crystal"):
                assert orientation.check_table(name, str(path)) is False


@pytest.mark.parametrize("n, wrong, missing", [(83, "orthogonal", "C16"), (88, "orthogonal", "C16"),
                                               (162, "trigonal_1", "C15")])
def test_known_gaps_are_reported_with_the_right_alternative(n, wrong, missing, tmp_path):
    path = tmp_path / "POSCAR"
    poscar.write(path, smallest_cell(n), "SG {}".format(n))
    with pytest.warns(UserWarning) as record:
        orientation.check_table(wrong, str(path))
    message = str(record[0].message)
    assert missing in message
    assert "{}_pre()".format(expected(n)[1]) in message and "auto_pre()" in message
