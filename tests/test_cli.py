"""Tests for check_symmetry(), the strains argument and the elastic-vasp command line."""
import contextlib
import io
import os
import warnings

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pytest  # noqa: E402

from elastic_vasp import cli, poscar  # noqa: E402
from elastic_vasp import elastic_gen as eg  # noqa: E402
from elastic_vasp.symmetry import check_symmetry, classify  # noqa: E402
from test_elastic import CASES, CUBIC, EV_PER_A3_TO_GPA, _poscar, outcar_text  # noqa: E402

EXAMPLES = os.path.join(os.path.dirname(__file__), os.pardir, "example")


def _run(argv):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        code = cli.main(argv)
    return code, out.getvalue(), err.getvalue()


# --- check_symmetry -------------------------------------------------------------

@pytest.mark.parametrize("folder, group, named", [
    ("Cubic", "Fm-3m (225)", "cubic"),
    ("Hexagonal", "P6_3/mmc (194)", "hexagonal"),
    ("Trigonal", "R-3c (167)", "trigonal_1"),
])
def test_check_symmetry_on_the_examples(folder, group, named):
    report = check_symmetry(os.path.join(EXAMPLES, folder, "POSCAR"), verbose=False)
    assert report.space_group == group
    assert report.recommended == named
    assert report.recommended_works  # the examples are in the standard orientation
    assert "works in this orientation" in report.summary()


def test_relaxation_noise_does_not_break_the_checks(tmp_path):
    """VASP-relaxed cells are symmetric only to ~1e-7 (as in example/Hexagonal)."""
    atoms = poscar.read(os.path.join(EXAMPLES, "Hexagonal", "POSCAR"))[0]
    rng = np.random.default_rng(0)
    atoms.set_cell(atoms.cell.array + rng.normal(0, 1e-7, (3, 3)), scale_atoms=True)
    path = tmp_path / "POSCAR"
    poscar.write(path, atoms, "noisy")
    report = check_symmetry(str(path), verbose=False)
    assert report.number == 194 and report.recommended_works
    from elastic_vasp import orientation
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert orientation.check_table("hexagonal", str(path))


def test_check_symmetry_reports_a_misoriented_crystal(tmp_path):
    from scipy.spatial.transform import Rotation
    atoms = poscar.read(os.path.join(EXAMPLES, "Hexagonal", "POSCAR"))[0]
    atoms.set_cell(atoms.cell.array @ Rotation.from_euler("y", 90, degrees=True).as_matrix().T,
                   scale_atoms=True)
    path = tmp_path / "POSCAR"
    poscar.write(path, atoms, "c along x")
    report = check_symmetry(str(path), verbose=False)
    assert report.recommended == "hexagonal" and not report.recommended_works
    # c along x is still symmetric about x, y and z, so the orthorhombic table works
    assert report.working == ["orthogonal", "monoclinic", "triclinic"]
    assert "NOT in this orientation" in report.summary()
    assert sorted(os.listdir(tmp_path)) == ["POSCAR"]  # read-only: nothing written


def test_classify():
    assert classify(1)[2] == "triclinic"
    assert classify(83) == ("tetragonal II", "4/m", "tetragonal_2")
    assert classify(162)[2] == "trigonal_2"   # '31m' type
    assert classify(166)[2] == "trigonal_1"
    with pytest.raises(ValueError):
        classify(231)


# --- strains passed directly ------------------------------------------------------

def test_symmetric_strains():
    assert eg.symmetric_strains(0.02, 5) == [-0.02, -0.01, 0.0, 0.01, 0.02]
    assert len(eg.symmetric_strains(0.03, 7)) == 7
    with pytest.raises(ValueError):
        eg.symmetric_strains(0.02, 2)


def test_strains_argument_replaces_strain_dat(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "POSCAR").write_text(_poscar(CUBIC))
    with contextlib.redirect_stdout(io.StringIO()):
        eg.cubic_pre(strains=eg.symmetric_strains(0.02, 5))
    assert sorted(os.listdir(tmp_path / "C44")) == [
        "strain_-0.01", "strain_-0.02", "strain_0.00", "strain_0.01", "strain_0.02"]
    assert not (tmp_path / "strain.dat").exists()  # not created when strains are given


@pytest.mark.parametrize("strains, match", [([0.01, 0.0, -0.0], "more than once"),
                                            ([-0.01, 0.01], "only 2 strain")])
def test_bad_strains_argument_touches_nothing(tmp_path, monkeypatch, strains, match):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "POSCAR").write_text(_poscar(CUBIC))
    with pytest.raises(ValueError, match=match):
        eg.cubic_pre(strains=strains)
    assert sorted(os.listdir(tmp_path)) == ["POSCAR"]


# --- command line --------------------------------------------------------------------

def test_cli_full_run(tmp_path):
    (tmp_path / "POSCAR").write_text(_poscar(CUBIC))
    code, out, err = _run(["-C", str(tmp_path), "check"])
    assert code == 0 and "Fm-3m" not in out and "Pm-3m (221)" in out and err == ""

    code, out, err = _run(["-C", str(tmp_path), "pre", "cubic", "--max-strain", "0.02", "--points", "5"])
    assert code == 0 and "elastic-vasp post cubic" in out
    assert (tmp_path / "C44" / "strain_0.02" / "POSCAR").exists()

    volume = abs(np.linalg.det(CUBIC))
    for path in tmp_path.glob("*/strain_*/POSCAR"):
        f = np.linalg.solve(CUBIC, poscar.lattice(path))
        e = np.array([f[0, 0] - 1, f[1, 1] - 1, f[2, 2] - 1, 2 * f[1, 2], 2 * f[0, 2], 2 * f[0, 1]])
        energy = volume * 0.5 * e @ CASES["cubic"] @ e / EV_PER_A3_TO_GPA
        (path.parent / "OUTCAR").write_text(outcar_text(volume, energy))
    code, out, err = _run(["-C", str(tmp_path), "post", "cubic", "--output", "C.txt"])
    assert code == 0, err
    np.testing.assert_allclose(np.loadtxt(tmp_path / "C.txt"), CASES["cubic"], atol=0.2)


def test_cli_auto_and_strains(tmp_path):
    (tmp_path / "POSCAR").write_text(_poscar(CUBIC))
    code, out, _ = _run(["-C", str(tmp_path), "pre", "auto", "--strains", "-0.01", "0", "0.01"])
    assert code == 0 and "3 independent elastic constants" in out
    assert (tmp_path / "elastic_vasp_auto.json").exists()


def test_cli_errors_exit_with_code_1(tmp_path):
    code, _, err = _run(["-C", str(tmp_path), "pre", "cubic"])
    assert code == 1 and err.startswith("error: POSCAR not found")
    (tmp_path / "POSCAR").write_text(_poscar(CUBIC))
    assert _run(["-C", str(tmp_path), "pre", "cubic", "--strains", "-0.01", "0", "0.01"])[0] == 0
    code, _, err = _run(["-C", str(tmp_path), "pre", "cubic", "--strains", "-0.01", "0", "0.01"])
    assert code == 1 and "already exist" in err
    code, _, err = _run(["-C", str(tmp_path), "post", "cubic"])
    assert code == 1 and "OUTCAR not found" in err


def test_cli_warnings_are_short(tmp_path):
    (tmp_path / "POSCAR").write_text(_poscar(np.array([[3.0, 0, 0], [0, 3.2, 0], [0.3, 0, 3.5]])))
    code, _, err = _run(["-C", str(tmp_path), "pre", "cubic", "--strains", "-0.01", "0", "0.01"])
    assert code == 0 and err.startswith("warning: cubic_pre() cannot measure")


def test_cli_systems_and_version():
    code, out, _ = _run(["systems"])
    assert code == 0 and "tetragonal_2" in out and "C11_C66_C16" in out
    with pytest.raises(SystemExit):
        _run(["--version"])


def test_cli_post_auto_standard_frame(tmp_path):
    from ase.lattice.cubic import FaceCenteredCubic
    from elastic_vasp import auto
    atoms = FaceCenteredCubic("Al", directions=[[1, -1, 0], [1, 1, -2], [1, 1, 1]],
                              latticeconstant=4.04, pbc=True)
    poscar.write(tmp_path / "POSCAR", atoms, "oriented")
    assert _run(["-C", str(tmp_path), "pre", "auto", "--max-strain", "0.02", "--points", "5"])[0] == 0
    dirs = np.array([[1, -1, 0], [1, 1, -2], [1, 1, 1]], float)
    m = auto.bond_matrix(dirs / np.linalg.norm(dirs, axis=1)[:, None])
    c_true = m @ CASES["cubic"] @ m.T
    ref = atoms.cell.array
    volume = abs(np.linalg.det(ref))
    for path in tmp_path.glob("*/strain_*/POSCAR"):
        f = np.linalg.solve(ref, poscar.lattice(path))
        e = np.array([f[0, 0] - 1, f[1, 1] - 1, f[2, 2] - 1, 2 * f[1, 2], 2 * f[0, 2], 2 * f[0, 1]])
        (path.parent / "OUTCAR").write_text(outcar_text(volume, volume * 0.5 * e @ c_true @ e / EV_PER_A3_TO_GPA))
    code, _, err = _run(["-C", str(tmp_path), "post", "auto", "--frame", "standard", "-o", "C.txt"])
    assert code == 0, err
    np.testing.assert_allclose(np.loadtxt(tmp_path / "C.txt"), CASES["cubic"], atol=0.3)
    code, _, err = _run(["-C", str(tmp_path), "post", "cubic", "--frame", "standard"])
    assert code == 1 and "only for 'auto'" in err
