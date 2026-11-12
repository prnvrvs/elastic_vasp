"""Tests for the crystal-symmetry orientation checks used by the *_pre functions."""
import contextlib
import io
import warnings

import numpy as np
import pytest

from elastic_vasp import elastic_gen as eg
from elastic_vasp import orientation

SQ3 = np.sqrt(3)
A, C = 3.0, 10.0

GOOD = {
    "monoclinic_b_unique": (orientation.check_monoclinic,
                            [[3.0, 0.0, 0.0], [0.0, 3.2, 0.0], [0.3, 0.0, 3.5]]),
    # base-centred monoclinic primitive cell, 2-fold still along y
    "monoclinic_c2m_primitive": (orientation.check_monoclinic,
                                 [[1.5, -1.6, 0.0], [1.5, 1.6, 0.0], [0.3, 0.0, 3.5]]),
    "hexagonal_c_along_z": (orientation.check_threefold_z,
                            [[A, 0, 0], [-A / 2, A * SQ3 / 2, 0], [0, 0, C]]),
    "rhombohedral_primitive_c_along_z": (orientation.check_threefold_z,
                                         [[A / 2, -A / (2 * SQ3), C / 3],
                                          [0, A / SQ3, C / 3],
                                          [-A / 2, -A / (2 * SQ3), C / 3]]),
}

BAD = {
    "monoclinic_c_unique": (orientation.check_monoclinic,
                            [[3.0, 0.0, 0.0], [0.4, 3.2, 0.0], [0.0, 0.0, 3.5]]),
    "hexagonal_c_along_x": (orientation.check_threefold_z,
                            [[0, A, 0], [0, -A / 2, A * SQ3 / 2], [C, 0, 0]]),
    "rhombohedral_111": (orientation.check_threefold_z,
                         [[3.0, 1.0, 1.0], [1.0, 3.0, 1.0], [1.0, 1.0, 3.0]]),
}


def _write_poscar(path, lattice, mode="Direct"):
    """One atom at the origin, so the crystal has the full symmetry of its lattice."""
    rows = "".join(" ".join("%.12f" % v for v in row) + "\n" for row in lattice)
    path.write_text("test\n1.0\n" + rows + "Al\n1\n" + mode + "\n0 0 0\n")


def _lattice(path):
    lines = path.read_text().splitlines()
    return np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)])


@pytest.mark.parametrize("name", sorted(GOOD))
def test_correct_orientation_passes_silently(name, tmp_path):
    check, lattice = GOOD[name]
    poscar = tmp_path / "POSCAR"
    _write_poscar(poscar, lattice)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert check(str(poscar)) is True
    assert not (tmp_path / orientation.REORIENTED_NAME).exists()


@pytest.mark.parametrize("name", sorted(BAD))
def test_misaligned_cell_warns_and_writes_rotated_copy(name, tmp_path):
    check, lattice = BAD[name]
    poscar = tmp_path / "POSCAR"
    _write_poscar(poscar, lattice)
    with pytest.warns(UserWarning, match="Continuing with the current orientation"):
        assert check(str(poscar)) is False

    reoriented = tmp_path / orientation.REORIENTED_NAME
    assert reoriented.exists()
    new = _lattice(reoriented)
    old = np.array(lattice, float)
    # pure rotation: same metric (lengths and angles), same volume sign
    np.testing.assert_allclose(new @ new.T, old @ old.T, atol=1e-9)
    assert np.linalg.det(new) == pytest.approx(np.linalg.det(old))
    # species, order and fractional coordinates are unchanged (only the cell rotates)
    from ase.io.vasp import read_vasp
    before, after = read_vasp(str(poscar)), read_vasp(str(reoriented))
    assert after.get_chemical_symbols() == before.get_chemical_symbols()
    np.testing.assert_allclose(after.get_scaled_positions(), before.get_scaled_positions(), atol=1e-12)
    # and the rotated copy passes the check
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert check(str(reoriented)) is True


def test_cartesian_poscar_also_gets_a_rotated_copy(tmp_path):
    check, lattice = BAD["monoclinic_c_unique"]
    poscar = tmp_path / "POSCAR"
    _write_poscar(poscar, lattice, mode="Cartesian")
    with pytest.warns(UserWarning, match="rotated copy"):
        check(str(poscar))
    reoriented = tmp_path / orientation.REORIENTED_NAME
    assert "Direct" in reoriented.read_text()
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert check(str(reoriented)) is True


def test_pre_warns_but_still_generates_folders(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _write_poscar(tmp_path / "POSCAR", BAD["rhombohedral_111"][1])
    (tmp_path / "strain.dat").write_text("-0.01,0.00,0.01\n")
    with pytest.warns(UserWarning, match="3-fold"):
        with contextlib.redirect_stdout(io.StringIO()):
            eg.trigonal_2_pre()
    assert (tmp_path / "C11_C44_C15" / "strain_0.01" / "POSCAR").exists()
    assert (tmp_path / orientation.REORIENTED_NAME).exists()


# --- crystal symmetry, not just the lattice ---------------------------------------

def _write_atoms(path, atoms):
    from elastic_vasp import poscar
    poscar.write(path, atoms, "test")


def _p2_crystal(cell, unique):
    """P2 crystal: atoms at general positions related only by a 2-fold about `unique`."""
    from ase import Atoms
    flip = {"x": (1, -1, -1), "y": (-1, 1, -1), "z": (-1, -1, 1)}[unique]
    general = [[0.1, 0.2, 0.3], [0.37, 0.61, 0.12]]
    frac = general + [list(np.multiply(p, flip)) for p in general]
    return Atoms("AlOAlO", cell=cell, pbc=True, scaled_positions=frac)


def _axes_along(atoms, order, target):
    return any(abs(abs(a @ np.asarray(target, float)) - 1) < 1e-6
               for a in orientation.symmetry_axes(atoms, order))


def test_pseudo_orthorhombic_cell_with_unique_axis_along_z_is_caught(tmp_path):
    # all angles 90 degrees: the lattice has 2-folds along x, y and z, the crystal only along z
    atoms = _p2_crystal(np.diag([3.0, 4.0, 5.0]), "z")
    src = tmp_path / "POSCAR"
    _write_atoms(src, atoms)
    with pytest.warns(UserWarning, match=r"space group P2 \(3\).*along Cartesian \(0.000, 0.000, 1.000\)"):
        assert orientation.check_monoclinic(str(src)) is False
    from elastic_vasp import poscar
    fixed = poscar.read(tmp_path / orientation.REORIENTED_NAME)[0]
    assert _axes_along(fixed, 2, [0, 1, 0])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert orientation.check_monoclinic(str(tmp_path / orientation.REORIENTED_NAME))


def test_monoclinic_crystal_with_unique_axis_along_y_passes(tmp_path):
    beta = np.radians(103)
    cell = [[3.0, 0, 0], [0, 4.0, 0], [5.0 * np.cos(beta), 0, 5.0 * np.sin(beta)]]
    src = tmp_path / "POSCAR"
    _write_atoms(src, _p2_crystal(cell, "y"))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert orientation.check_monoclinic(str(src)) is True


def test_structure_without_the_axis_says_so(tmp_path):
    from ase import Atoms
    triclinic = Atoms("AlO", cell=[[3, 0.2, 0.1], [0.3, 4, 0.2], [0.4, 0.1, 5]], pbc=True,
                      scaled_positions=[[0, 0, 0], [0.31, 0.17, 0.43]])
    src = tmp_path / "POSCAR"
    _write_atoms(src, triclinic)
    with pytest.warns(UserWarning, match="no 2-fold \\(unique\\) axis at all"):
        assert orientation.check_monoclinic(str(src)) is False
    assert not (tmp_path / orientation.REORIENTED_NAME).exists()


def test_cubic_crystal_in_trigonal_check_reports_higher_symmetry(tmp_path):
    src = tmp_path / "POSCAR"
    _write_poscar(src, [[0, 2, 2], [2, 0, 2], [2, 2, 0]])  # fcc primitive cell: 3-folds along <111>
    with pytest.warns(UserWarning, match="It has 4 such axes"):
        assert orientation.check_threefold_z(str(src), system="trigonal") is False


@pytest.mark.parametrize("seed", range(5))
def test_rotated_corundum_is_put_back_with_its_3fold_axis_along_z(tmp_path, seed):
    import os
    from scipy.spatial.transform import Rotation
    from elastic_vasp import poscar
    example = os.path.join(os.path.dirname(__file__), os.pardir, "example", "Trigonal", "POSCAR")
    atoms = poscar.read(example)[0]
    group = orientation.space_group(atoms)
    rot = Rotation.random(random_state=seed).as_matrix()
    atoms.set_cell(atoms.cell.array @ rot.T, scale_atoms=True)
    src = tmp_path / "POSCAR"
    _write_atoms(src, atoms)
    with pytest.warns(UserWarning):
        assert orientation.check_threefold_z(str(src), system="trigonal") is False
    fixed = poscar.read(tmp_path / orientation.REORIENTED_NAME)[0]
    assert orientation.space_group(fixed) == group == "R-3c (167)"
    assert _axes_along(fixed, 3, [0, 0, 1])
    np.testing.assert_allclose(fixed.get_scaled_positions(), atoms.get_scaled_positions(), atol=1e-9)
    assert np.linalg.det(fixed.cell.array) > 0
