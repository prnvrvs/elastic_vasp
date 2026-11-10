"""Regression tests for the energy-strain pipeline.

Each test generates the deformation folders with the real ``*_pre`` function,
writes synthetic OUTCARs whose energy is E = V * 1/2 * eps^T C eps for a known
stiffness matrix C, runs the matching ``*_post`` function and checks that C is
recovered.
"""
import contextlib
import io
import warnings

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pytest  # noqa: E402

from elastic_vasp import coefficient, constants  # noqa: E402
from elastic_vasp import elastic_gen as eg  # noqa: E402

EV_PER_A3_TO_GPA = 160.2176621
# b-unique monoclinic cell (2-fold along y); fine for the systems without a check
GENERIC = np.array([[3.0, 0.0, 0.0], [0.0, 3.2, 0.0], [0.3, 0.0, 3.5]])
HEXAGONAL = np.array([[3.0, 0.0, 0.0], [-1.5, 1.5 * np.sqrt(3), 0.0], [0.0, 0.0, 5.0]])


def _poscar(lattice):
    return "test\n1.0\n" + "".join(
        " ".join("%.12f" % v for v in row) + "\n" for row in lattice
    ) + "Al\n1\nDirect\n0 0 0\n"


def _sym(upper):
    m = np.array(upper, float)
    return np.triu(m) + np.triu(m, 1).T


def _voigt_strain(poscar_path, reference):
    """Recover the applied Voigt strain from a generated POSCAR."""
    with open(poscar_path) as fh:
        lines = fh.readlines()
    lattice = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)])
    f = np.linalg.solve(reference, lattice)  # lattice = reference @ f
    return np.array([f[0, 0] - 1, f[1, 1] - 1, f[2, 2] - 1,
                     2 * f[1, 2], 2 * f[0, 2], 2 * f[0, 1]])


def _run(system, c_true, tmp_path, monkeypatch):
    lattice = HEXAGONAL if system in ("hexagonal", "trigonal_1", "trigonal_2") else GENERIC
    volume = abs(np.linalg.det(lattice))
    monkeypatch.chdir(tmp_path)
    (tmp_path / "POSCAR").write_text(_poscar(lattice))
    (tmp_path / "strain.dat").write_text("-0.02,-0.01,0.00,0.01,0.02\n")
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # correctly oriented cells must not warn
        with contextlib.redirect_stdout(io.StringIO()):
            getattr(eg, system + "_pre")()

    for poscar in tmp_path.glob("*/strain_*/POSCAR"):
        eps = _voigt_strain(poscar, lattice)
        energy = volume * 0.5 * eps @ c_true @ eps / EV_PER_A3_TO_GPA
        (poscar.parent / "OUTCAR").write_text(
            "  volume of cell :   %.8f\n"
            "  free  energy   TOTEN  =   %.12f eV\n" % (volume, energy)
        )

    with contextlib.redirect_stdout(io.StringIO()):
        c_out = getattr(constants, system + "_post")()
    return c_out


def _c66(c11, c12):
    return (c11 - c12) / 2


CASES = {
    "cubic": _sym([[250, 120, 120, 0, 0, 0], [0, 250, 120, 0, 0, 0], [0, 0, 250, 0, 0, 0],
                   [0, 0, 0, 80, 0, 0], [0, 0, 0, 0, 80, 0], [0, 0, 0, 0, 0, 80]]),
    "hexagonal": _sym([[300, 100, 80, 0, 0, 0], [0, 300, 80, 0, 0, 0], [0, 0, 350, 0, 0, 0],
                       [0, 0, 0, 90, 0, 0], [0, 0, 0, 0, 90, 0], [0, 0, 0, 0, 0, _c66(300, 100)]]),
    "orthogonal": _sym([[300, 110, 90, 0, 0, 0], [0, 280, 95, 0, 0, 0], [0, 0, 260, 0, 0, 0],
                        [0, 0, 0, 70, 0, 0], [0, 0, 0, 0, 75, 0], [0, 0, 0, 0, 0, 85]]),
    "trigonal_1": _sym([[300, 100, 80, 20, 0, 0], [0, 300, 80, -20, 0, 0], [0, 0, 350, 0, 0, 0],
                        [0, 0, 0, 90, 0, 0], [0, 0, 0, 0, 90, 20], [0, 0, 0, 0, 0, _c66(300, 100)]]),
    # trigonal II adds C15, C25 = -C15, C46 = -C15
    "trigonal_2": _sym([[300, 100, 80, 20, 15, 0], [0, 300, 80, -20, -15, 0], [0, 0, 350, 0, 0, 0],
                        [0, 0, 0, 90, 0, -15], [0, 0, 0, 0, 90, 20], [0, 0, 0, 0, 0, _c66(300, 100)]]),
    "monoclinic": _sym([[300, 110, 90, 0, 12, 0], [0, 280, 95, 0, 8, 0], [0, 0, 260, 0, 6, 0],
                        [0, 0, 0, 70, 0, 10], [0, 0, 0, 0, 75, 0], [0, 0, 0, 0, 0, 85]]),
    "triclinic": (lambda a: a @ a.T / 6 + np.eye(6) * 150)(
        np.random.default_rng(0).normal(size=(6, 6)) * 30),
}


@pytest.mark.parametrize("system", sorted(CASES))
def test_recovers_stiffness(system, tmp_path, monkeypatch):
    c_true = CASES[system]
    c_out = _run(system, c_true, tmp_path, monkeypatch)
    np.testing.assert_allclose(c_out, c_true, atol=0.2)


def test_pre_refuses_to_overwrite_existing_results(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "POSCAR").write_text(_poscar(GENERIC))
    (tmp_path / "strain.dat").write_text("-0.01,0.00,0.01\n")
    eg.cubic_pre()
    outcar = tmp_path / "C44" / "strain_0.01" / "OUTCAR"
    outcar.write_text("finished run")

    with pytest.raises(FileExistsError):
        eg.cubic_pre()
    assert outcar.exists()

    eg.cubic_pre(overwrite=True)
    assert not outcar.exists()


def test_oszicar_prefers_e0(tmp_path):
    oszicar = tmp_path / "OSZICAR"
    oszicar.write_text(
        "DAV:   1    -0.1E+02\n"
        "   1 F= -.10000000E+02 E0= -.11000000E+02  d E =-.1E+02\n"
    )
    assert coefficient._read_last_oszicar_energy(str(oszicar)) == pytest.approx(-11.0)


def test_coeff_ignores_non_strain_folders(tmp_path):
    for label, energy in (("-0.01", 1e-3), ("0.00", 0.0), ("0.01", 1e-3)):
        folder = tmp_path / ("strain_" + label)
        folder.mkdir()
        (folder / "OUTCAR").write_text(
            "  volume of cell :   10.0\n  free  energy   TOTEN  =   %.6f eV\n" % energy
        )
    (tmp_path / "notes").mkdir()
    # E/V = a * s^2 with a = 1e-3 / 10 / 1e-4 eV/A^3
    assert coefficient.coeff(str(tmp_path)) == pytest.approx(1.0 * EV_PER_A3_TO_GPA, rel=1e-6)
