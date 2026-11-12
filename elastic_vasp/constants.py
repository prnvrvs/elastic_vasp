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
from elastic_vasp import auto
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


def analyse(system, base_dir=None, after_matrix=None):
    """Post-process the deformation folders of `system` in `base_dir` (default: cwd).

    `system` is a name from systems.SYSTEMS or a systems.CrystalSystem.
    """
    spec = SYSTEMS[system] if isinstance(system, str) else system
    base_dir = os.getcwd() if base_dir is None else base_dir
    with outcar_cache():  # each OUTCAR is scanned once for the whole run
        return _analyse(spec, base_dir, after_matrix)


def _analyse(spec, base_dir, after_matrix=None):
    _check_inputs(spec, base_dir)
    art.loading()

    coefficients = {d.folder: _coeff_from(base_dir, d.folder) for d in spec.deformations}

    print("second order coefficients (a in E/V0 = a*delta^2, GPa) are as follow:")
    for d in spec.deformations:
        print("{}/= {}".format(d.folder, coefficients[d.folder]))

    print("\n    ===== Methodology =====")
    print("    strain (e1..e6) per deformation, and a = 1/2 e^T C e:")
    design = spec.design_matrix()
    for k, d in enumerate(spec.deformations):
        print("    {:<16} {:<54} a = {}".format(
            d.folder, d.describe(), spec.energy_expression(k, design)))
    print("    solved as a linear system for: " + ", ".join(spec.constants))

    values, stiffness = spec.solve(coefficients)
    print()
    for name in spec.constants:
        print("{}= {}".format(name, values[name]))

    print('''
    ==== stiffness matrix ====''')
    print(np.around(stiffness, decimals=2))
    if after_matrix is not None:
        after_matrix(stiffness)

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
    """Cubic: elastic constants from the finished VASP runs.

    Fits the VASP energies in the cubic_pre() folders, solves for C11, C12, C44,
    prints the stiffness matrix, the stability check and the Voigt-Reuss-Hill
    averages, and returns the 6x6 stiffness matrix in GPa.
    """
    return analyse("cubic")


def tetragonal_1_post():
    """Tetragonal I: elastic constants from the finished VASP runs.

    Fits the VASP energies in the tetragonal_1_pre() folders, solves for C11, C12,
    C13, C33, C44, C66, prints the stiffness matrix, the stability check and the
    Voigt-Reuss-Hill averages, and returns the 6x6 stiffness matrix in GPa.
    """
    return analyse("tetragonal_1")


def tetragonal_2_post():
    """Tetragonal II: elastic constants from the finished VASP runs.

    Fits the VASP energies in the tetragonal_2_pre() folders, solves for C11, C12,
    C13, C33, C44, C66, C16, prints the stiffness matrix, the stability check and
    the Voigt-Reuss-Hill averages, and returns the 6x6 stiffness matrix in GPa.
    """
    return analyse("tetragonal_2")


def hexagonal_post():
    """Hexagonal: elastic constants from the finished VASP runs.

    Fits the VASP energies in the hexagonal_pre() folders, solves for C11, C12, C13,
    C33, C44, prints the stiffness matrix, the stability check and the Voigt-Reuss-
    Hill averages, and returns the 6x6 stiffness matrix in GPa.
    """
    return analyse("hexagonal")


def orthogonal_post():
    """Orthorhombic: elastic constants from the finished VASP runs.

    Fits the VASP energies in the orthogonal_pre() folders, solves for C11, C22,
    C33, C12, C13, C23, C44, C55, C66, prints the stiffness matrix, the stability
    check and the Voigt-Reuss-Hill averages, and returns the 6x6 stiffness matrix in
    GPa.
    """
    return analyse("orthogonal")


def trigonal_1_post():
    """Trigonal I: elastic constants from the finished VASP runs.

    Fits the VASP energies in the trigonal_1_pre() folders, solves for C11, C12,
    C13, C33, C44, C14, prints the stiffness matrix, the stability check and the
    Voigt-Reuss-Hill averages, and returns the 6x6 stiffness matrix in GPa.
    """
    return analyse("trigonal_1")


def trigonal_2_post():
    """Trigonal II: elastic constants from the finished VASP runs.

    Fits the VASP energies in the trigonal_2_pre() folders, solves for C11, C12,
    C13, C33, C44, C14, C15, prints the stiffness matrix, the stability check and
    the Voigt-Reuss-Hill averages, and returns the 6x6 stiffness matrix in GPa.
    """
    return analyse("trigonal_2")


def monoclinic_post():
    """Monoclinic: elastic constants from the finished VASP runs.

    Fits the VASP energies in the monoclinic_pre() folders, solves for C11, C22,
    C33, C12, C13, C23, C44, C55, C66, C15, C25, C35, C46, prints the stiffness
    matrix, the stability check and the Voigt-Reuss-Hill averages, and returns the
    6x6 stiffness matrix in GPa.
    """
    return analyse("monoclinic")


def triclinic_post():
    """Triclinic: elastic constants from the finished VASP runs.

    Fits the VASP energies in the triclinic_pre() folders, solves for C11, C12, C13,
    C14, C15, C16, C22, C23, C24, C25, C26, C33, C34, C35, C36, C44, C45, C46, C55,
    C56, C66, prints the stiffness matrix, the stability check and the Voigt-Reuss-
    Hill averages, and returns the 6x6 stiffness matrix in GPa.
    """
    return analyse("triclinic")


def auto_post(base_dir=None, frame="input"):
    """Post-process folders made by elastic_gen.auto_pre(), using its saved setup.

    Prints the stiffness matrix in the orientation of the input POSCAR and, in
    addition, in the crystal's standard setting (cube axes along x, y, z;
    hexagonal c along z; ...) together with the textbook constants of its crystal
    system (e.g. C11, C12, C44 for cubic). Both describe the same crystal; the
    stability check and the polycrystalline averages are the same in either frame.

    Returns the 6x6 matrix in GPa: in the input orientation (``frame="input"``,
    default) or in the standard setting (``frame="standard"``).
    """
    if frame not in ("input", "standard"):
        raise ValueError("frame must be 'input' or 'standard', got {!r}".format(frame))
    base_dir = os.getcwd() if base_dir is None else base_dir
    path = os.path.join(base_dir, auto.SETUP_FILE)
    if not os.path.exists(path):
        raise FileNotFoundError("{} not found: run elastic_gen.auto_pre() first (it saves the "
                                "deformations to use).".format(path))
    standard = auto.load_standard(path)
    if standard is None and frame == "standard":
        raise ValueError("{} was written by an older version without the standard setting; "
                         "rerun auto_pre(overwrite=True)".format(path))
    results = {}

    def report(c_input):  # printed right after the stiffness matrix
        if standard is not None:
            results["standard"] = _report_standard(c_input, *standard)

    c_input = analyse(auto.load(path), base_dir, after_matrix=report)
    return results["standard"] if frame == "standard" else c_input


def _report_standard(c_input, rotation, number):
    from elastic_vasp.symmetry import classify

    c_standard = auto.to_standard(c_input, rotation)
    system, laue, _ = classify(number)
    name, values, residual = auto.standard_constants(c_standard, number)
    print("""
    ==== in the crystal's standard setting ====
    space group {} - {} (Laue class {}); rotation from your cell to the standard frame:""".format(
        number, system, laue))
    print(np.array2string(rotation, precision=4, suppress_small=True, prefix="    "))
    print("    stiffness matrix (GPa):")
    print(np.around(c_standard, decimals=2))
    print("    {} constants (as from {}_post()):".format(system, name))
    for key, value in values.items():
        print("{}= {:.4f}".format(key, value))
    if residual > 0.5:
        print(_RED + " warning: the matrix deviates by {:.2f} GPa from the {} form; check the "
              "symmetry tolerance and the calculations".format(residual, system) + _RESET)
    return c_standard
