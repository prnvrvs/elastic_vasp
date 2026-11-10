"""Tests for the lattice orientation checks used by the *_pre functions."""
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
    rows = "".join(" ".join("%.12f" % v for v in row) + "\n" for row in lattice)
    path.write_text("test\n1.0\n" + rows + "Al Si\n1 1\n" + mode + "\n0 0 0\n0.25 0.5 0.1\n")


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
    # everything after the lattice (species, counts, fractional coordinates) is unchanged
    assert reoriented.read_text().splitlines()[5:] == poscar.read_text().splitlines()[5:]
    # and the rotated copy passes the check
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert check(str(reoriented)) is True


def test_cartesian_poscar_warns_without_writing_copy(tmp_path):
    check, lattice = BAD["monoclinic_c_unique"]
    poscar = tmp_path / "POSCAR"
    _write_poscar(poscar, lattice, mode="Cartesian")
    with pytest.warns(UserWarning, match="Cartesian coordinates"):
        check(str(poscar))
    assert not (tmp_path / orientation.REORIENTED_NAME).exists()


def test_pre_warns_but_still_generates_folders(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _write_poscar(tmp_path / "POSCAR", BAD["rhombohedral_111"][1])
    (tmp_path / "strain.dat").write_text("-0.01,0.00,0.01\n")
    with pytest.warns(UserWarning, match="3-fold"):
        with contextlib.redirect_stdout(io.StringIO()):
            eg.trigonal_2_pre()
    assert (tmp_path / "C11_C44_C15" / "strain_0.01" / "POSCAR").exists()
    assert (tmp_path / orientation.REORIENTED_NAME).exists()
