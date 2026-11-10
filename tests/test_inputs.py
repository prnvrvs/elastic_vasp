"""Tests for POSCAR scale handling, volume lookup, strain.dat parsing and fit checks."""
import contextlib
import io
import os
import warnings

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pytest  # noqa: E402

from elastic_vasp import coefficient, orientation, poscar  # noqa: E402
from elastic_vasp import elastic_gen as eg  # noqa: E402
from elastic_vasp.target import calc_gen  # noqa: E402

CUBE = np.diag([3.0, 3.0, 3.0])


def _write(path, scale_line, rows, mode="Direct"):
    body = "".join(" ".join("%.10f" % v for v in row) + "\n" for row in rows)
    path.write_text("test\n" + scale_line + "\n" + body + "Al\n1\n" + mode + "\n0 0 0\n")


def _outcar(volume, e0):
    return ("  volume of cell :   %.8f\n"
            "  energy  without entropy=   %.12f  energy(sigma->0) =   %.12f\n"
            " General timing and accounting informations for this job:\n" % (volume, e0, e0))


# --- POSCAR scale factors ----------------------------------------------------

@pytest.mark.parametrize("scale_line, rows, actual", [
    ("1.0", CUBE, CUBE),
    ("2.0", np.eye(3) * 1.5, CUBE),
    ("-27.0", np.eye(3) * 7.0, CUBE),             # volume 27 -> 3 A cube
    ("3.0 1.5 0.5", np.diag([1.0, 2.0, 6.0]), CUBE),  # VASP 6 per-axis scales
    ("  -27.0   ! target volume", np.eye(3), CUBE),  # trailing comment
])
def test_actual_lattice(tmp_path, scale_line, rows, actual):
    path = tmp_path / "POSCAR"
    _write(path, scale_line, rows)
    np.testing.assert_allclose(poscar.lattice(poscar.read(path)), actual)
    assert coefficient._read_volume_from_poscar(str(path)) == pytest.approx(27.0)


@pytest.mark.parametrize("scale_line, rows", [("-27.0", np.eye(3) * 7.0),
                                              ("3.0 1.5 0.5", np.diag([1.0, 2.0, 6.0]))])
def test_strained_poscar_keeps_the_strain_with_special_scales(tmp_path, scale_line, rows):
    src, out = tmp_path / "POSCAR", tmp_path / "out"
    _write(src, scale_line, rows)
    s = 0.02
    calc_gen(s, s, s, 0, 0, 0, poscar_path=str(src), output_path=str(out))
    lines = poscar.read(out)
    assert lines[1].strip() == "1.0"  # VASP would otherwise rescale and undo the strain
    np.testing.assert_allclose(poscar.lattice(lines), CUBE * (1 + s), atol=1e-10)
    assert lines[5:] == poscar.read(src)[5:]


def test_plain_scale_is_kept_as_written(tmp_path):
    src, out = tmp_path / "POSCAR", tmp_path / "out"
    _write(src, "2.0", np.eye(3) * 1.5)
    calc_gen(0.01, 0, 0, 0, 0, 0, poscar_path=str(src), output_path=str(out))
    lines = poscar.read(out)
    assert lines[1] == "2.0\n"
    np.testing.assert_allclose(poscar.lattice(lines), np.diag([3.03, 3.0, 3.0]))


def test_orientation_copy_with_negative_scale(tmp_path):
    a, c = 3.0, 5.0
    hex_c_along_x = [[0, a, 0], [0, -a / 2, a * np.sqrt(3) / 2], [c, 0, 0]]
    volume = abs(np.linalg.det(hex_c_along_x))
    src = tmp_path / "POSCAR"
    _write(src, "-%.10f" % volume, np.array(hex_c_along_x) * 2.0)  # rows twice too long
    with pytest.warns(UserWarning):
        orientation.check_threefold_z(str(src))
    reoriented = poscar.read(tmp_path / orientation.REORIENTED_NAME)
    assert reoriented[1].strip() == "1.0"
    lat = poscar.lattice(reoriented)
    assert abs(np.linalg.det(lat)) == pytest.approx(volume)
    assert orientation.check_threefold_z(str(tmp_path / orientation.REORIENTED_NAME))


# --- reference volume -----------------------------------------------------------

def _strain_folders(base, labels, volume=27.0):
    for label in labels:
        folder = base / ("strain_" + label)
        folder.mkdir(parents=True)
        (folder / "OUTCAR").write_text(_outcar(volume, 1e-3 * float(label) ** 2 * 1e4))


def test_coeff_accepts_relative_base_dir(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _strain_folders(tmp_path / "C44", ["-0.01", "0.00", "0.01"])
    assert coefficient.coeff("C44") == pytest.approx(coefficient.coeff(str(tmp_path / "C44")))


def test_coeff_without_zero_strain_uses_unstrained_poscar(tmp_path):
    _write(tmp_path / "POSCAR", "1.0", CUBE)
    _strain_folders(tmp_path / "C44", ["-0.02", "-0.01", "0.01", "0.02"], volume=99.0)
    value = coefficient.coeff(str(tmp_path / "C44"))
    # E = 1e-3 eV at s = 0.01, V0 = 27 A^3 from ../POSCAR (not 99 from the strained runs)
    assert value == pytest.approx(1e-3 / 27.0 / 1e-4 * 160.2176621, rel=1e-6)


def test_coeff_without_any_reference_volume_says_so(tmp_path):
    _strain_folders(tmp_path / "C44", ["-0.02", "-0.01", "0.01"])
    with pytest.raises(FileNotFoundError, match="reference volume"):
        coefficient.coeff(str(tmp_path / "C44"))


def test_flat_energies_warn(tmp_path):
    for label in ("-0.01", "0.00", "0.01"):
        folder = tmp_path / ("strain_" + label)
        folder.mkdir()
        (folder / "OUTCAR").write_text(_outcar(27.0, -5.0))
    with pytest.warns(UserWarning, match="identical"):
        coefficient.coeff(str(tmp_path))


# --- strain.dat -------------------------------------------------------------------

def test_strain_dat_reads_every_line(tmp_path):
    path = tmp_path / "strain.dat"
    path.write_text("# strains\n-0.02, -0.01\n0.00 0.01  # zero and +1%\n\n0.02\n")
    assert eg._read_strains(str(path)) == ["-0.02", "-0.01", "0.00", "0.01", "0.02"]


@pytest.mark.parametrize("text, match", [
    ("0.01, 0.00, -0.00\n", "more than once: -0.00"),
    ("0.01, 0.0l, 0.02\n", "'0.0l' is not a number"),
    ("-0.01, 0.01\n", "only 2 strain value"),
    ("# nothing\n\n", "does not contain"),
])
def test_bad_strain_dat_fails_before_touching_folders(tmp_path, monkeypatch, text, match):
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / "POSCAR", "1.0", CUBE)
    (tmp_path / "strain.dat").write_text(text)
    with pytest.raises(ValueError, match=match):
        with contextlib.redirect_stdout(io.StringIO()):
            eg.cubic_pre()
    assert not os.path.exists(tmp_path / "C44")


def test_missing_strain_dat_is_created_with_defaults(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / "POSCAR", "1.0", CUBE)
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        eg.cubic_pre()
    assert (tmp_path / "strain.dat").read_text() == eg.DEFAULT_STRAINS + "\n"
    assert "strain.dat not found: created it" in out.getvalue()
    assert sorted(os.listdir(tmp_path / "C44")) == sorted(
        "strain_" + s for s in eg.DEFAULT_STRAINS.split(","))


# --- unfinished runs, too few points, POSCAR checks, OUTCAR reads ---------------

def test_killed_run_is_reported(tmp_path):
    folder = tmp_path / "strain_0.01"
    folder.mkdir()
    (folder / "OUTCAR").write_text(
        "  volume of cell :   27.0\n"
        "  energy  without entropy=   -5.0  energy(sigma->0) =   -5.0\n")  # no footer
    with pytest.raises(coefficient.MissingEnergyError, match="Run did not finish"):
        coefficient._read_energy(str(folder))


def test_coeff_refuses_fewer_than_three_points(tmp_path):
    _write(tmp_path / "POSCAR", "1.0", CUBE)
    _strain_folders(tmp_path / "C44", ["-0.01", "0.01"])
    with pytest.raises(ValueError, match="needs at least 3"):
        coefficient.coeff(str(tmp_path / "C44"))


@pytest.mark.parametrize("pre", ["cubic_pre", "hexagonal_pre", "monoclinic_pre", "trigonal_2_pre"])
def test_missing_poscar_touches_nothing_even_with_overwrite(tmp_path, monkeypatch, pre):
    monkeypatch.chdir(tmp_path)
    _write(tmp_path / "POSCAR", "1.0", CUBE)
    (tmp_path / "strain.dat").write_text("-0.01,0.00,0.01\n")
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")  # orientation warnings for the cube are irrelevant here
        getattr(eg, pre)()
    before = sorted(str(p.relative_to(tmp_path)) for p in tmp_path.rglob("*"))
    (tmp_path / "POSCAR").unlink()
    with pytest.raises(FileNotFoundError, match="nothing was changed"):
        with contextlib.redirect_stdout(io.StringIO()):
            getattr(eg, pre)(overwrite=True)
    after = sorted(str(p.relative_to(tmp_path)) for p in tmp_path.rglob("*"))
    assert after == [p for p in before if p != "POSCAR"]


def test_unreadable_poscar_touches_nothing(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    (tmp_path / "POSCAR").write_text("test\nnot-a-number\n")
    (tmp_path / "strain.dat").write_text("-0.01,0.00,0.01\n")
    with pytest.raises(ValueError, match="could not be read"):
        eg.cubic_pre()
    assert not (tmp_path / "C44").exists()


def test_each_outcar_is_read_once_during_post(tmp_path, monkeypatch):
    import builtins
    from elastic_vasp import constants

    monkeypatch.chdir(tmp_path)
    _write(tmp_path / "POSCAR", "1.0", CUBE)
    (tmp_path / "strain.dat").write_text("-0.01,0.00,0.01\n")
    with contextlib.redirect_stdout(io.StringIO()):
        eg.cubic_pre()
    outcars = list(tmp_path.glob("*/strain_*"))
    for folder in outcars:
        s = float(folder.name.split("_", 1)[1])
        (folder / "OUTCAR").write_text(_outcar(27.0, 1e-3 * (s * 100) ** 2))

    opened = {}
    real_open = builtins.open

    def counting_open(file, *args, **kwargs):
        if str(file).endswith("OUTCAR"):
            opened[str(file)] = opened.get(str(file), 0) + 1
        return real_open(file, *args, **kwargs)

    monkeypatch.setattr(builtins, "open", counting_open)
    with contextlib.redirect_stdout(io.StringIO()):
        constants.cubic_post()
    monkeypatch.setattr(builtins, "open", real_open)
    assert len(opened) == len(outcars)
    assert set(opened.values()) == {1}


def _relax_outcar(nsw, ibrion, converged):
    text = ("   NSW    =     %d    number of steps for IOM\n"
            "   IBRION =      %d    ionic relax: 0-MD 1-quasi-New 2-CG\n"
            "  volume of cell :   27.0\n"
            "  energy  without entropy=   -5.0  energy(sigma->0) =   -5.0\n" % (nsw, ibrion))
    if converged:
        text += " reached required accuracy - stopping structural energy minimisation\n"
    return text + " General timing and accounting informations for this job:\n"


@pytest.mark.parametrize("nsw, ibrion, converged, ok", [
    (40, 2, True, True),     # relaxation converged
    (40, 2, False, False),   # ran out of NSW steps
    (40, 1, False, False),
    (0, 2, False, True),     # static run: nothing to converge
    (40, -1, False, True),   # IBRION = -1: no ionic update
    (40, 0, False, True),    # MD is not a relaxation
])
def test_relaxation_must_converge(tmp_path, nsw, ibrion, converged, ok):
    folder = tmp_path / "strain_0.01"
    folder.mkdir()
    (folder / "OUTCAR").write_text(_relax_outcar(nsw, ibrion, converged))
    if ok:
        assert coefficient._read_energy(str(folder)) == pytest.approx(-5.0)
    else:
        with pytest.raises(coefficient.MissingEnergyError, match="relaxation did not converge"):
            coefficient._read_energy(str(folder))


def test_trigonal_example_is_correctly_oriented():
    import warnings as _w
    example = os.path.join(os.path.dirname(__file__), os.pardir, "example", "Trigonal", "POSCAR")
    with _w.catch_warnings():
        _w.simplefilter("error")
        assert orientation.check_threefold_z(example, system="trigonal")
