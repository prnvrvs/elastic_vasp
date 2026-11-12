"""Golden-output test: pre/post must reproduce the recorded reference exactly.

tests/data/reference.json was generated from the hand-written per-system code
(before the move to the table-driven implementation in systems.py). It holds,
for every crystal system,

* the lattice and fractional atom positions of every generated POSCAR (so the
  deformations written to disk are unchanged; compared numerically so that
  ASE's text formatting, which can vary between versions, does not matter), and
* the 6x6 stiffness matrix returned by *_post() for synthetic energies.

Regenerate (only when an output change is intended) with:

    python tests/test_reference.py
"""
import contextlib
import io
import json
import os
import sys
import tempfile
import warnings
from pathlib import Path

import numpy as np
import pytest
from ase.io.vasp import read_vasp

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from test_elastic import (CASES, EV_PER_A3_TO_GPA, GENERIC, HEXAGONAL, TETRAGONAL,  # noqa: E402
                          _poscar, _voigt_strain, outcar_text)

from elastic_vasp import constants  # noqa: E402
from elastic_vasp import elastic_gen as eg  # noqa: E402

REFERENCE = Path(__file__).resolve().parent / "data" / "reference.json"
STRAINS = "-0.02,-0.01,0.00,0.01,0.02\n"


@contextlib.contextmanager
def _chdir(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def generate(system, workdir):
    workdir = Path(workdir)
    # the lattices the reference was recorded with (kept as they were; GENERIC is
    # monoclinic, so the symmetry check warns for cubic/orthogonal - ignored below)
    if system in ("hexagonal", "trigonal_1", "trigonal_2"):
        lattice = HEXAGONAL
    elif system.startswith("tetragonal"):
        lattice = TETRAGONAL
    else:
        lattice = GENERIC
    volume = abs(np.linalg.det(lattice))
    c_true = CASES[system]
    with _chdir(workdir), warnings.catch_warnings(), contextlib.redirect_stdout(io.StringIO()):
        warnings.simplefilter("error")
        warnings.filterwarnings("ignore", message=".*cannot measure this crystal")
        (workdir / "POSCAR").write_text(_poscar(lattice))
        (workdir / "strain.dat").write_text(STRAINS)
        getattr(eg, system + "_pre")()

        structures = {}
        for poscar in sorted(workdir.glob("*/strain_*/POSCAR")):
            atoms = read_vasp(str(poscar))
            structures[str(poscar.relative_to(workdir))] = {
                "cell": atoms.cell.array.tolist(),
                "fractional": atoms.get_scaled_positions(wrap=False).tolist(),
            }
            eps = _voigt_strain(poscar, lattice)
            energy = volume * 0.5 * eps @ c_true @ eps / EV_PER_A3_TO_GPA
            (poscar.parent / "OUTCAR").write_text(outcar_text(volume, energy))
        c_out = getattr(constants, system + "_post")()
    return {"structures": structures, "C": np.asarray(c_out).tolist()}


@pytest.mark.parametrize("system", sorted(CASES))
def test_matches_reference(system, tmp_path):
    reference = json.loads(REFERENCE.read_text())[system]
    result = generate(system, tmp_path)
    assert sorted(result["structures"]) == sorted(reference["structures"])
    for path, ref in reference["structures"].items():
        got = result["structures"][path]
        np.testing.assert_allclose(got["cell"], ref["cell"], rtol=0, atol=1e-9, err_msg=path)
        np.testing.assert_allclose(got["fractional"], ref["fractional"], rtol=0, atol=1e-9, err_msg=path)
    np.testing.assert_allclose(result["C"], reference["C"], rtol=0, atol=1e-8)


if __name__ == "__main__":
    data = {}
    for name in sorted(CASES):
        with tempfile.TemporaryDirectory() as tmp:
            data[name] = generate(name, tmp)
    REFERENCE.parent.mkdir(exist_ok=True)
    REFERENCE.write_text(json.dumps(data, indent=1, sort_keys=True) + "\n")
    print("wrote", REFERENCE, {k: len(v["structures"]) for k, v in data.items()})
