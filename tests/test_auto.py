"""Tests for the orientation-aware setup (auto.py, auto_pre, auto_post)."""
import contextlib
import io
import os
import warnings

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pytest  # noqa: E402
from ase import Atoms  # noqa: E402
from scipy.spatial.transform import Rotation  # noqa: E402

from elastic_vasp import auto, constants, poscar  # noqa: E402
from elastic_vasp import elastic_gen as eg  # noqa: E402
from elastic_vasp.systems import SYSTEMS  # noqa: E402
from test_elastic import CASES, EV_PER_A3_TO_GPA, HEXAGONAL, outcar_text  # noqa: E402

EXAMPLE = os.path.join(os.path.dirname(__file__), os.pardir, "example", "Trigonal", "POSCAR")


def _single(cell):
    return Atoms("Al", cell=cell, pbc=True, scaled_positions=[[0, 0, 0]])


def _corundum():
    return poscar.read(EXAMPLE)[0]


def _ilmenite():
    """Corundum with Fe/Ti ordered along c, keeping the inversion centre: R-3."""
    atoms = _corundum()
    symbols = atoms.get_chemical_symbols()
    symbols[0] = symbols[2] = "Fe"
    symbols[1] = symbols[3] = "Ti"
    atoms.set_chemical_symbols(symbols)
    return atoms


def _p2(cell, unique="y"):
    flip = {"y": (-1, 1, -1), "z": (-1, -1, 1)}[unique]
    general = [[0.1, 0.2, 0.3], [0.37, 0.61, 0.12]]
    frac = general + [list(np.multiply(p, flip)) for p in general]
    return Atoms("AlOAlO", cell=cell, pbc=True, scaled_positions=frac)


BETA = np.radians(103)
MONO_CELL = [[3.0, 0, 0], [0, 4.0, 0], [5.0 * np.cos(BETA), 0, 5.0 * np.sin(BETA)]]
STANDARD = {
    "cubic": _single(np.eye(3) * 3),
    "hexagonal": _single(HEXAGONAL),
    "orthogonal": _single(np.diag([3.0, 4.0, 5.0])),
    "trigonal_1": _corundum(),
    "trigonal_2": _ilmenite(),
    "monoclinic": _p2(MONO_CELL),
    "triclinic": Atoms("AlO", cell=[[3, 0.2, 0.1], [0.3, 4, 0.2], [0.4, 0.1, 5]], pbc=True,
                       scaled_positions=[[0, 0, 0], [0.31, 0.17, 0.43]]),
}


@pytest.mark.parametrize("name", sorted(STANDARD))
def test_standard_orientation_gives_the_textbook_constants(name):
    spec, _ = auto.from_structure(STANDARD[name])
    table = SYSTEMS[name]
    assert sorted(spec.constants) == sorted(table.constants)
    assert len(spec.deformations) == len(table.constants)
    # same relations between entries (signs, C66 = (C11-C12)/2, ...)
    values = dict(zip(table.constants, np.random.default_rng(0).uniform(10, 90, len(table.constants))))
    np.testing.assert_allclose(spec.matrix([values[n] for n in spec.constants]),
                               table.matrix([values[n] for n in table.constants]), atol=1e-9)


@pytest.mark.parametrize("name", sorted(STANDARD))
def test_constants_are_invariant_under_the_crystal_symmetry(name):
    atoms = STANDARD[name].copy()
    atoms.set_cell(atoms.cell.array @ Rotation.random(random_state=3).as_matrix().T, scale_atoms=True)
    spec, basis = auto.from_structure(atoms)
    assert len(spec.constants) == len(SYSTEMS[name].constants)  # same count in any orientation
    for r in auto.crystal_rotations(atoms):
        m = auto.bond_matrix(r)
        for b in basis.values():
            np.testing.assert_allclose(m @ b @ m.T, b, atol=1e-9)


def _run_auto(atoms, c_true, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    poscar.write(tmp_path / "POSCAR", atoms, "test")
    (tmp_path / "strain.dat").write_text("-0.02,-0.01,0.00,0.01,0.02\n")
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("error")
        eg.auto_pre()
    reference = atoms.cell.array
    volume = abs(np.linalg.det(reference))
    for path in tmp_path.glob("*/strain_*/POSCAR"):
        f = np.linalg.solve(reference, poscar.lattice(path))
        eps = np.array([f[0, 0] - 1, f[1, 1] - 1, f[2, 2] - 1, 2 * f[1, 2], 2 * f[0, 2], 2 * f[0, 1]])
        energy = volume * 0.5 * eps @ c_true @ eps / EV_PER_A3_TO_GPA
        (path.parent / "OUTCAR").write_text(outcar_text(volume, energy))
    with contextlib.redirect_stdout(io.StringIO()):
        return constants.auto_post()


@pytest.mark.parametrize("name, seed", [("trigonal_1", 1), ("trigonal_2", 2), ("hexagonal", 3),
                                        ("monoclinic", 4), ("cubic", 5)])
def test_rotated_crystal_end_to_end(name, seed, tmp_path, monkeypatch):
    """Crystal rotated arbitrarily: auto_pre/auto_post recover C in the rotated frame."""
    atoms = STANDARD[name].copy()
    rot = Rotation.random(random_state=seed).as_matrix()
    atoms.set_cell(atoms.cell.array @ rot.T, scale_atoms=True)
    m = auto.bond_matrix(rot)
    c_true = m @ CASES[name] @ m.T
    c_out = _run_auto(atoms, c_true, tmp_path, monkeypatch)
    np.testing.assert_allclose(c_out, c_true, atol=0.3)


def test_pseudo_orthorhombic_monoclinic_end_to_end(tmp_path, monkeypatch):
    """All cell angles 90 degrees, unique axis along z: no rotation needed with auto."""
    atoms = _p2(np.diag([3.0, 4.0, 5.0]), unique="z")
    rot = Rotation.from_euler("x", 90, degrees=True).as_matrix()  # y -> z
    m = auto.bond_matrix(rot)
    c_true = m @ CASES["monoclinic"] @ m.T
    c_out = _run_auto(atoms, c_true, tmp_path, monkeypatch)
    np.testing.assert_allclose(c_out, c_true, atol=0.3)
    assert abs(c_out[0, 5]) > 1  # C16 appears in this frame (it is C15 in the standard one)


def test_auto_post_needs_the_saved_setup(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError, match="run elastic_gen.auto_pre"):
        constants.auto_post()


def test_refused_auto_pre_keeps_the_old_setup(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    poscar.write(tmp_path / "POSCAR", STANDARD["cubic"], "test")
    (tmp_path / "strain.dat").write_text("-0.01,0.00,0.01\n")
    with contextlib.redirect_stdout(io.StringIO()):
        eg.auto_pre()
    saved = (tmp_path / auto.SETUP_FILE).read_text()
    poscar.write(tmp_path / "POSCAR", STANDARD["hexagonal"], "other structure")
    with contextlib.redirect_stdout(io.StringIO()), pytest.raises(FileExistsError):
        eg.auto_pre()
    assert (tmp_path / auto.SETUP_FILE).read_text() == saved


def test_setup_round_trip(tmp_path):
    spec, basis = auto.from_structure(_ilmenite())
    auto.save(tmp_path / "setup.json", spec, basis, auto.SYMPREC)
    loaded = auto.load(tmp_path / "setup.json")
    assert loaded.constants == spec.constants
    assert loaded.deformations == spec.deformations
    values = np.arange(1, len(spec.constants) + 1, dtype=float)
    np.testing.assert_allclose(loaded.matrix(values), spec.matrix(values))


@pytest.mark.parametrize("name", sorted(STANDARD))
def test_deformations_are_well_conditioned_in_any_orientation(name):
    """Errors in the fitted coefficients must not be blown up in C (tables: 1.9 - 5.1)."""
    for seed in range(4):
        atoms = STANDARD[name].copy()
        atoms.set_cell(atoms.cell.array @ Rotation.random(random_state=seed).as_matrix().T,
                       scale_atoms=True)
        assert auto.amplification(auto.from_structure(atoms)[0]) < 8


def test_noisy_energies_give_bounded_errors(tmp_path, monkeypatch):
    """Rotated corundum with 0.5 meV noise on every energy: C stays close to the truth."""
    atoms = STANDARD["trigonal_1"].copy()
    rot = Rotation.random(random_state=11).as_matrix()
    atoms.set_cell(atoms.cell.array @ rot.T, scale_atoms=True)
    m = auto.bond_matrix(rot)
    c_true = m @ CASES["trigonal_1"] @ m.T
    rng = np.random.default_rng(5)
    real_outcar_text = outcar_text

    def noisy(volume, energy):
        return real_outcar_text(volume, energy + rng.normal(0, 5e-4))
    monkeypatch.setitem(globals(), "outcar_text", noisy)
    c_out = _run_auto(atoms, c_true, tmp_path, monkeypatch)
    assert np.abs(c_out - c_true).max() < 15  # GPa; with 80000x amplification this was hopeless


def test_standard_frame_printout_shows_the_pattern(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    poscar.write(tmp_path / "POSCAR", STANDARD["hexagonal"], "hex")
    (tmp_path / "strain.dat").write_text("-0.01,0.00,0.01\n")
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        eg.auto_pre()
    text = out.getvalue()
    assert "P6/mmm (191)" in text and "5 independent elastic constants" in text
    assert "0.5*C11 - 0.5*C12" in text  # C66 shown in terms of C11, C12
    assert "Error amplification" in text


def test_oriented_fcc_supercell_x110_y112_z111(tmp_path, monkeypatch):
    """Common MD-style orientation: 6-atom fcc cell with x || [1-10], y || [11-2], z || [111]."""
    from ase.lattice.cubic import FaceCenteredCubic
    from elastic_vasp.symmetry import check_symmetry

    atoms = FaceCenteredCubic("Al", directions=[[1, -1, 0], [1, 1, -2], [1, 1, 1]],
                              latticeconstant=4.04, pbc=True)
    spec, _ = auto.from_structure(atoms)
    assert len(spec.constants) == 3 and len(spec.deformations) == 3
    path = tmp_path / "POSCAR"
    poscar.write(path, atoms, "fcc 110/112/111")
    report = check_symmetry(str(path), verbose=False)
    assert report.n_constants == 3 and not report.recommended_works
    assert report.working[:2] == ["trigonal_1", "trigonal_2"]  # 3-fold along z, 2-fold along x

    dirs = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]], float)
    m = auto.bond_matrix(dirs / np.linalg.norm(dirs, axis=1)[:, None])
    c_true = m @ CASES["cubic"] @ m.T
    np.testing.assert_allclose(_run_auto(atoms, c_true, tmp_path, monkeypatch), c_true, atol=0.3)


# --- standard setting ---------------------------------------------------------------

def _run_auto_frames(atoms, c_true, tmp_path, monkeypatch):
    c_input = _run_auto(atoms, c_true, tmp_path, monkeypatch)
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        c_standard = constants.auto_post(frame="standard")
    return c_input, c_standard, out.getvalue()


@pytest.mark.parametrize("name, seed", [("cubic", 21), ("hexagonal", 22), ("trigonal_1", 23),
                                        ("trigonal_2", 24), ("monoclinic", 25)])
def test_standard_setting_after_arbitrary_rotation(name, seed, tmp_path, monkeypatch):
    atoms = STANDARD[name].copy()
    rot = Rotation.random(random_state=seed).as_matrix()
    atoms.set_cell(atoms.cell.array @ rot.T, scale_atoms=True)
    m = auto.bond_matrix(rot)
    c_true = m @ CASES[name] @ m.T
    c_input, c_standard, text = _run_auto_frames(atoms, c_true, tmp_path, monkeypatch)
    np.testing.assert_allclose(c_input, c_true, atol=0.3)
    assert "in the crystal's standard setting" in text
    # same crystal: frame-independent eigenvalues
    np.testing.assert_allclose(np.linalg.eigvalsh(c_standard), np.linalg.eigvalsh(CASES[name]), atol=0.3)
    # in the textbook form of its crystal system
    number = orientation_number(atoms)
    table_name, values, residual = auto.standard_constants(c_standard, number)
    assert residual < 0.3
    if name in ("cubic", "hexagonal"):     # no choice of in-plane axes: values are unique
        np.testing.assert_allclose(c_standard, CASES[name], atol=0.3)
    if name.startswith("trigonal"):        # C14 / C15 signs depend on which equivalent a is chosen
        for key in ("C11", "C12", "C13", "C33", "C44"):
            assert values[key] == pytest.approx(CASES[name][int(key[1]) - 1, int(key[2]) - 1], abs=0.3)


def orientation_number(atoms):
    return auto.standard_setting(atoms)[1]


def test_oriented_fcc_gives_textbook_cubic_constants(tmp_path, monkeypatch):
    from ase.lattice.cubic import FaceCenteredCubic
    atoms = FaceCenteredCubic("Al", directions=[[1, -1, 0], [1, 1, -2], [1, 1, 1]],
                              latticeconstant=4.04, pbc=True)
    dirs = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]], float)
    m = auto.bond_matrix(dirs / np.linalg.norm(dirs, axis=1)[:, None])
    _, c_standard, text = _run_auto_frames(atoms, m @ CASES["cubic"] @ m.T, tmp_path, monkeypatch)
    np.testing.assert_allclose(c_standard, CASES["cubic"], atol=0.3)  # 250 / 120 / 80
    assert "cubic constants (as from cubic_post())" in text and "C44= 80.0" in text


def test_old_setup_without_standard_setting(tmp_path, monkeypatch):
    import json
    atoms = STANDARD["cubic"]
    c_true = CASES["cubic"]
    _run_auto(atoms, c_true, tmp_path, monkeypatch)
    data = json.loads((tmp_path / auto.SETUP_FILE).read_text())
    data.pop("standard_rotation"); data.pop("space_group_number")
    (tmp_path / auto.SETUP_FILE).write_text(json.dumps(data))
    with contextlib.redirect_stdout(io.StringIO()):
        np.testing.assert_allclose(constants.auto_post(), c_true, atol=0.3)
        with pytest.raises(ValueError, match="older version"):
            constants.auto_post(frame="standard")
