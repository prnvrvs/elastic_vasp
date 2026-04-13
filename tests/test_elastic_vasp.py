import os
import shutil
import tempfile
import unittest
from pathlib import Path

import numpy as np

from elastic_vasp import constants as c
from elastic_vasp import elastic_gen as eg
from elastic_vasp.coefficient import coeff


REPO_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_POSCAR = REPO_ROOT / "example" / "Cubic" / "POSCAR"
VOL = 100.0
ENERGY_FACTOR = 160.2176621


def _write_fake_vasp_outputs(base_dir, folders, coeff_value=2.0):
    for folder in folders:
        folder_dir = base_dir / folder
        for strain_dir in sorted(folder_dir.iterdir()):
            if not strain_dir.is_dir():
                continue
            strain = float(strain_dir.name.split("_", 1)[1])
            energy = coeff_value * (strain ** 2) * VOL / ENERGY_FACTOR
            (strain_dir / "OSZICAR").write_text(
                " 1 F= {:.12f} E0= {:.12f} d E = {:.12f}\n".format(energy, energy, energy)
            )
            (strain_dir / "OUTCAR").write_text("volume of cell is {:.6f}\n".format(VOL))


def _write_fake_outcar_energy(strain_dir, energy):
    strain_dir.mkdir(parents=True, exist_ok=True)
    strain_dir.joinpath("OUTCAR").write_text(
        " free  energy   TOTEN  = {: .8f} eV\n".format(energy)
    )


def _write_poscar_with_volume(path, volume):
    lattice = np.array([[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, volume / 100.0]])
    lines = [
        "test\n",
        "1.0\n",
        "{:16.10f} {:16.10f} {:16.10f}\n".format(*lattice[0]),
        "{:16.10f} {:16.10f} {:16.10f}\n".format(*lattice[1]),
        "{:16.10f} {:16.10f} {:16.10f}\n".format(*lattice[2]),
        "Al\n",
        "1\n",
        "Direct\n",
        "0 0 0\n",
    ]
    path.write_text("".join(lines))


class ElasticVaspWorkflowTests(unittest.TestCase):
    def _run_in_tempdir(self):
        workdir = Path(tempfile.mkdtemp(prefix="elastic-vasp-test-"))
        self.addCleanup(lambda: shutil.rmtree(workdir, ignore_errors=True))
        return workdir

    def test_preprocessing_creates_expected_folders(self):
        cases = [
            (eg.cubic_pre, ["C11_C12_I", "C11_C12_II", "C44"]),
            (eg.orthogonal_pre, ["C11_C22_C12_I", "C11_C22_C12_II", "C22_C33_C23_I", "C22_C33_C23_II", "C11_C33_C13_I", "C11_C33_C13_II", "C44", "C55", "C66"]),
            (eg.hexagonal_pre, ["C11_C12_I", "C11_C12_II", "C11_C33_C13_I", "C44", "C11_C12_C13_C33"]),
            (eg.trigonal_1_pre, ["C11_C12_I", "C11_C12_II", "C11_C33_C13_I", "C11_C33_C13_II", "C11_C44_C14_I", "C11_C44_C14_II"]),
            (eg.trigonal_2_pre, ["C11_C12_I", "C11_C12_II", "C11_C33_C13_I", "C11_C33_C13_II", "C11_C44_C14_I", "C11_C44_C14_II", "C11_C44_C15"]),
            (eg.monoclinic_pre, ["C11_C22_C12_I", "C11_C22_C12_II", "C22_C33_C23_I", "C22_C33_C23_II", "C11_C33_C13_I", "C11_C33_C13_II", "C44", "C11_C55_C15_I", "C11_C55_C15_II", "C22_C55_C25", "C33_C55_C35", "C44_C66_C46_I", "C44_C66_C46_II"]),
            (eg.triclinic_pre, ["C11", "C11_C22_C12", "C22", "C22_C33_C23", "C33", "C11_C33_C13", "C44", "C55", "C66", "C11_C44_C14", "C11_C55_C15", "C11_C66_C16", "C22_C44_C24", "C22_C55_C25", "C22_C66_C26", "C33_C44_C34", "C33_C55_C35", "C33_C66_C36", "C44_C55_C45", "C44_C66_C46", "C55_C66_C56"]),
        ]

        for pre_fn, folders in cases:
            workdir = self._run_in_tempdir()
            shutil.copy(EXAMPLE_POSCAR, workdir / "POSCAR")
            (workdir / "strain.dat").write_text("-0.01,0.000,0.01\n")
            cwd = os.getcwd()
            try:
                os.chdir(workdir)
                pre_fn()
                self.assertEqual(sorted(p.name for p in workdir.iterdir() if p.is_dir()), sorted(folders))
                for folder in folders:
                    strain_dirs = sorted((workdir / folder).iterdir())
                    self.assertEqual(len(strain_dirs), 3)
                    self.assertTrue((workdir / folder / "strain_0.00" / "POSCAR").exists())
            finally:
                os.chdir(cwd)

    def test_coeff_detects_zero_folder(self):
        workdir = self._run_in_tempdir()
        cwd = os.getcwd()
        try:
            os.chdir(workdir)
            for label, energy in [("-0.01", 1.0), ("0.000", 0.5), ("0.01", 1.0)]:
                strain_dir = workdir / f"strain_{label}"
                strain_dir.mkdir()
                (strain_dir / "OSZICAR").write_text(
                    " 1 F= {:.6f} E0= {:.6f} d E = {:.6f}\n".format(energy, energy, energy)
                )
                (strain_dir / "OUTCAR").write_text("volume of cell is 100.000\n")
            value = coeff(base_dir=str(workdir))
            self.assertTrue(np.isfinite(value))
            self.assertGreater(value, 0.0)
        finally:
            os.chdir(cwd)

    def test_coeff_prefers_outcar_and_poscar_volume(self):
        workdir = self._run_in_tempdir()
        cwd = os.getcwd()
        try:
            os.chdir(workdir)
            for label in ["-0.01", "0.000", "0.01"]:
                strain_dir = workdir / f"strain_{label}"
                strain_dir.mkdir()
                strain = float(label)
                energy = 2.5 * (strain ** 2) * VOL / ENERGY_FACTOR
                _write_fake_outcar_energy(strain_dir, energy)
            _write_poscar_with_volume(workdir / "strain_0.000" / "POSCAR", VOL)
            value = coeff(base_dir=str(workdir))
            self.assertTrue(np.isfinite(value))
            self.assertGreater(value, 0.0)
            self.assertAlmostEqual(float((workdir / "Result.txt").read_text().splitlines()[1].split()[1]), 0.0, places=7)
        finally:
            os.chdir(cwd)

    def test_all_postprocessing_functions_run(self):
        cases = [
            (eg.cubic_pre, c.cubic_post, ["C11_C12_I", "C11_C12_II", "C44"]),
            (eg.orthogonal_pre, c.orthogonal_post, ["C11_C22_C12_I", "C11_C22_C12_II", "C22_C33_C23_I", "C22_C33_C23_II", "C11_C33_C13_I", "C11_C33_C13_II", "C44", "C55", "C66"]),
            (eg.hexagonal_pre, c.hexagonal_post, ["C11_C12_I", "C11_C12_II", "C11_C33_C13_I", "C44", "C11_C12_C13_C33"]),
            (eg.trigonal_1_pre, c.trigonal_1_post, ["C11_C12_I", "C11_C12_II", "C11_C33_C13_I", "C11_C33_C13_II", "C11_C44_C14_I", "C11_C44_C14_II"]),
            (eg.trigonal_2_pre, c.trigonal_2_post, ["C11_C12_I", "C11_C12_II", "C11_C33_C13_I", "C11_C33_C13_II", "C11_C44_C14_I", "C11_C44_C14_II", "C11_C44_C15"]),
            (eg.monoclinic_pre, c.monoclinic_post, ["C11_C22_C12_I", "C11_C22_C12_II", "C22_C33_C23_I", "C22_C33_C23_II", "C11_C33_C13_I", "C11_C33_C13_II", "C44", "C11_C55_C15_I", "C11_C55_C15_II", "C22_C55_C25", "C33_C55_C35", "C44_C66_C46_I", "C44_C66_C46_II"]),
            (eg.triclinic_pre, c.triclinic_post, ["C11", "C11_C22_C12", "C22", "C22_C33_C23", "C33", "C11_C33_C13", "C44", "C55", "C66", "C11_C44_C14", "C11_C55_C15", "C11_C66_C16", "C22_C44_C24", "C22_C55_C25", "C22_C66_C26", "C33_C44_C34", "C33_C55_C35", "C33_C66_C36", "C44_C55_C45", "C44_C66_C46", "C55_C66_C56"]),
        ]

        for pre_fn, post_fn, folders in cases:
            workdir = self._run_in_tempdir()
            shutil.copy(EXAMPLE_POSCAR, workdir / "POSCAR")
            (workdir / "strain.dat").write_text("-0.03,-0.02,0.000,0.02,0.03\n")
            cwd = os.getcwd()
            try:
                os.chdir(workdir)
                pre_fn()
                _write_fake_vasp_outputs(workdir, folders, coeff_value=2.0)
                mat = post_fn()
                self.assertEqual(mat.shape, (6, 6))
                self.assertTrue(np.all(np.isfinite(mat)))
            finally:
                os.chdir(cwd)


if __name__ == "__main__":
    unittest.main()
