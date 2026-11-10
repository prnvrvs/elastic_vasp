"""
Created on Sat Jun  5 19:16:26 2021

@author: pranav
"""
import os
import re
import shutil

import numpy as np

from elastic_vasp import orientation, poscar
from elastic_vasp.systems import SYSTEMS
from elastic_vasp.target import calc_gen


DEFAULT_STRAINS = "-0.03,-0.02,-0.01,0.00,0.01,0.02,0.03"


def _read_strains(path="strain.dat"):
    """All strain values in strain.dat: comma- and/or whitespace-separated, on
    any number of lines; '#' starts a comment.

    If the file does not exist it is created with DEFAULT_STRAINS (the values
    used in the examples), so the choice is visible and can be edited.
    """
    if not os.path.exists(path):
        with open(path, "w") as f:
            f.write(DEFAULT_STRAINS + "\n")
        print("{} not found: created it with the default strains {}. Edit it and "
              "rerun with overwrite=True to use other values.".format(path, DEFAULT_STRAINS))
    strains = []
    with open(path) as f:
        for line in f:
            line = line.split("#", 1)[0]
            strains.extend(item for item in re.split(r"[,\s]+", line) if item)
    if not strains:
        raise ValueError("strain.dat does not contain any strain values")
    # fail early (before any folder is touched) on typos and repeated values,
    # e.g. "0.00" and "-0.00" would both become strain_0.00
    seen = set()
    duplicates = []
    for item in strains:
        try:
            value = float(item)
        except ValueError:
            raise ValueError("{}: '{}' is not a number".format(path, item)) from None
        if value in seen:
            duplicates.append(item)
        seen.add(value)
    if duplicates:
        raise ValueError("strain.dat lists the same strain more than once: {}".format(
            ", ".join(duplicates)))
    if len(strains) < 3:
        raise ValueError("{} has only {} strain value(s); the quadratic E(s) fit needs at "
                         "least 3 (e.g. {})".format(path, len(strains), DEFAULT_STRAINS))
    return strains


def _normalize_strain_label(raw_strain):
    strain_value = float(raw_strain)
    if abs(strain_value) < 1e-12:
        return "0.00"
    return raw_strain.strip()


def _reset_folder(folder):
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    elif os.path.exists(folder):
        os.remove(folder)



def _check_poscar(path="POSCAR"):
    """Fail before any folder is touched if POSCAR is missing or unreadable."""
    if not os.path.exists(path):
        raise FileNotFoundError("{} not found in {}; nothing was changed.".format(path, os.getcwd()))
    try:
        cell = poscar.lattice(poscar.read(path))
    except (ValueError, IndexError) as exc:
        raise ValueError("{} could not be read ({}); nothing was changed.".format(path, exc)) from None
    if abs(np.linalg.det(cell)) < 1e-8:
        raise ValueError("{} has a degenerate lattice (zero volume); nothing was changed.".format(path))


def _generate_preprocessing(folder_specs, overwrite=False):
    _check_poscar()
    strains = _read_strains()
    cwd = os.getcwd()

    # Refuse to wipe existing deformation folders (they may hold finished
    # VASP runs) unless the caller explicitly asks for it.
    existing = [folder for folder, _ in folder_specs if os.path.exists(folder)]
    if existing and not overwrite:
        raise FileExistsError(
            "Deformation folders already exist: {}. They may contain finished "
            "VASP calculations. Remove them yourself or call this function "
            "with overwrite=True to regenerate them.".format(", ".join(existing))
        )

    for folder, deformation in folder_specs:
        _reset_folder(folder)
        os.mkdir(folder)

        for raw_strain in strains:
            strain_value = float(raw_strain)
            strain_label = _normalize_strain_label(raw_strain)
            strain_dir = os.path.join(folder, "strain_" + strain_label)

            os.mkdir(strain_dir)
            calc_gen(*deformation(strain_value))

            output_path = os.path.join(cwd, "output.txt")
            if not os.path.exists(output_path):
                raise IOError("calc_gen did not create output.txt")
            shutil.copy(output_path, os.path.join(strain_dir, "POSCAR"))
            os.remove(output_path)


def generate(system, overwrite=False):
    """Create the deformation folders for a crystal system from systems.SYSTEMS."""
    spec = SYSTEMS[system]
    _generate_preprocessing([(d.folder, d.at) for d in spec.deformations], overwrite)


def cubic_pre(overwrite=False):
    generate("cubic", overwrite)


def orthogonal_pre(overwrite=False):
    generate("orthogonal", overwrite)


def hexagonal_pre(overwrite=False):
    _check_poscar()
    orientation.check_threefold_z(system="hexagonal")
    generate("hexagonal", overwrite)


def trigonal_1_pre(overwrite=False):
    print("Note: trigonal_1 (point groups 32, 3m, -3m) also assumes x lies along a 2-fold "
          "axis / normal to a mirror plane. This holds for e.g. P-3m1, R-3m, P321, R3m with "
          "a along x, but NOT for '31m'-type groups (P-31m, P-31c, P31m, P31c, P312, ...). "
          "If unsure, use trigonal_2_pre() instead: it fits C14 and C15 together and works "
          "for any in-plane orientation.")
    _check_poscar()
    orientation.check_threefold_z(system="trigonal")
    generate("trigonal_1", overwrite)


def trigonal_2_pre(overwrite=False):
    _check_poscar()
    orientation.check_threefold_z(system="trigonal")
    generate("trigonal_2", overwrite)


def monoclinic_pre(overwrite=False):
    _check_poscar()
    orientation.check_monoclinic()
    generate("monoclinic", overwrite)


def triclinic_pre(overwrite=False):
    generate("triclinic", overwrite)
