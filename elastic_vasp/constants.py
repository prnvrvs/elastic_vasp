"""
Created on Sat Jun  5 19:16:26 2021

@author: pranav

Post-processing: fit E(s) in every deformation folder, solve for the
independent elastic constants of the crystal system (see systems.py), print
the results and return the 6x6 stiffness matrix in GPa (Voigt notation).
"""
import os
from glob import glob

import numpy as np

from elastic_vasp import art
from elastic_vasp import average
from elastic_vasp.coefficient import MissingEnergyError, check_energies, coeff, outcar_cache
from elastic_vasp.systems import SYSTEMS

_GREEN = "\033[1m\033[1;32m"
_RED = "\033[1m\033[1;31m"
_RESET = "\033[0m"


def _coeff_from(base_dir, folder):
    return coeff(os.path.join(base_dir, folder))


def _check_inputs(spec, base_dir):
    """Stop before any fitting if a deformation folder or an OUTCAR energy is missing."""
    problems = []
    strain_folders = []
    for d in spec.deformations:
        folder = os.path.join(base_dir, d.folder)
        found = sorted(glob(os.path.join(folder, "strain_*/")))
        if not found:
            problems.append("no strain_* folders in {}".format(folder))
        elif len(found) < 3:
            problems.append("only {} strain_* folder(s) in {}; the quadratic fit needs at "
                            "least 3".format(len(found), folder))
        strain_folders.extend(f.rstrip("/") for f in found)
    try:
        check_energies(strain_folders)
    except MissingEnergyError as exc:
        problems.extend(line.strip() for line in str(exc).splitlines()[1:])
    if problems:
        raise MissingEnergyError(
            "Cannot post-process {}: {} problem(s) found:\n  ".format(spec.name, len(problems))
            + "\n  ".join(problems))


def analyse(system, base_dir=None):
    """Post-process the deformation folders of `system` in `base_dir` (default: cwd)."""
    spec = SYSTEMS[system]
    base_dir = os.getcwd() if base_dir is None else base_dir
    with outcar_cache():  # each OUTCAR is scanned once for the whole run
        return _analyse(spec, base_dir)


def _analyse(spec, base_dir):
    _check_inputs(spec, base_dir)
    art.loading()

    coefficients = {d.folder: _coeff_from(base_dir, d.folder) for d in spec.deformations}

    print("second order coefficients (a in E/V = a*s^2, GPa) are as follow:")
    for d in spec.deformations:
        print("{}/= {}".format(d.folder, coefficients[d.folder]))

    print("\n    ===== Methodology =====")
    print("    strain (e1..e6) per deformation, and a = 1/2 e^T C e:")
    design = spec.design_matrix()
    for k, d in enumerate(spec.deformations):
        print("    {:<16} {:<44} a = {}".format(
            d.folder, d.describe(), spec.energy_expression(k, design)))
    print("    solved as a linear system for: " + ", ".join(spec.constants))

    values, stiffness = spec.solve(coefficients)
    print()
    for name in spec.constants:
        print("{}= {}".format(name, values[name]))

    print('''
    ==== stiffness matrix ====''')
    print(np.around(stiffness, decimals=2))

    print('''
    check stability criteria
    {}
         '''.format(spec.stability_text))
    if np.all(np.real(np.linalg.eigvals(stiffness)) > 0):
        print(_GREEN + " structure is mechanically stable" + _RESET)
    else:
        print(_RED + " structure is mechanically unstable" + _RESET)
    average.avg(stiffness)
    art.baginga()
    return stiffness


def cubic_post():
    return analyse("cubic")


def hexagonal_post():
    return analyse("hexagonal")


def orthogonal_post():
    return analyse("orthogonal")


def trigonal_1_post():
    return analyse("trigonal_1")


def trigonal_2_post():
    return analyse("trigonal_2")


def monoclinic_post():
    return analyse("monoclinic")


def triclinic_post():
    return analyse("triclinic")
