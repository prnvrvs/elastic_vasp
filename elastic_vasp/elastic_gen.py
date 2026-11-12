"""
Created on Sat Jun  5 19:16:26 2021

@author: pranav
"""
import os
import re
import shutil

import numpy as np

from elastic_vasp import auto, orientation, poscar
from elastic_vasp.systems import SYSTEMS
from elastic_vasp.symmetry import check_symmetry  # noqa: F401  (re-exported)
from elastic_vasp.target import calc_gen


DEFAULT_STRAINS = "-0.03,-0.02,-0.01,0.00,0.01,0.02,0.03"


def symmetric_strains(max_strain=0.03, points=7):
    """Evenly spaced strains from -max_strain to +max_strain (includes 0 for odd points).

    >>> symmetric_strains(0.02, 5)
    [-0.02, -0.01, 0.0, 0.01, 0.02]
    """
    if points < 3:
        raise ValueError("points must be at least 3 for the quadratic fit, got {}".format(points))
    if not 0 < max_strain < 0.5:
        raise ValueError("max_strain must be between 0 and 0.5, got {}".format(max_strain))
    return [round(float(v), 10) for v in np.linspace(-max_strain, max_strain, points)]


def _format_strain(value):
    """Folder label for a strain given as a number: 0.01 -> '0.01', 0 -> '0.00'."""
    if abs(value) < 1e-12:
        return "0.00"
    return "{:.10g}".format(value)


def _validate_strains(strains, source):
    """Check strain labels before any folder is touched; returns them unchanged."""
    if not strains:
        raise ValueError("{} does not contain any strain values".format(source))
    # typos and repeated values (e.g. "0.00" and "-0.00" would both become strain_0.00)
    seen = set()
    duplicates = []
    for item in strains:
        try:
            value = float(item)
        except ValueError:
            raise ValueError("{}: '{}' is not a number".format(source, item)) from None
        if value in seen:
            duplicates.append(item)
        seen.add(value)
    if duplicates:
        raise ValueError("{} lists the same strain more than once: {}".format(
            source, ", ".join(duplicates)))
    if len(strains) < 3:
        raise ValueError("{} has only {} strain value(s); the quadratic E(s) fit needs at "
                         "least 3 (e.g. {})".format(source, len(strains), DEFAULT_STRAINS))
    return strains


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
    return _validate_strains(strains, path)


def _strains_from_argument(strains):
    """Strain labels from a list of numbers (or numeric strings) passed by the caller."""
    labels = []
    for value in strains:
        if isinstance(value, str):
            labels.append(value.strip())
        else:
            labels.append(_format_strain(float(value)))
    return _validate_strains(labels, "strains argument")


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
        cell = poscar.lattice(path)
    except ValueError as exc:
        raise ValueError("{}; nothing was changed.".format(exc)) from None
    if abs(np.linalg.det(cell)) < 1e-8:
        raise ValueError("{} has a degenerate lattice (zero volume); nothing was changed.".format(path))


def _check_table(system):
    """Warn if the named system's deformations cannot measure this crystal (see orientation.check_table)."""
    return orientation.check_table(system)


def _generate_preprocessing(folder_specs, overwrite=False, strains=None):
    _check_poscar()
    strains = _read_strains() if strains is None else _strains_from_argument(strains)
    poscar_path = os.path.abspath("POSCAR")

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
            # write straight into the strain folder (no temporary file in the cwd)
            calc_gen(*deformation(strain_value), poscar_path=poscar_path,
                     output_path=os.path.join(strain_dir, "POSCAR"))


def generate(system, overwrite=False, strains=None):
    """Create the deformation folders for a crystal system.

    `system` is a name from systems.SYSTEMS or a systems.CrystalSystem.
    """
    spec = SYSTEMS[system] if isinstance(system, str) else system
    _generate_preprocessing([(d.folder, d.at) for d in spec.deformations], overwrite, strains)


def auto_pre(overwrite=False, symprec=auto.SYMPREC, strains=None):
    """Choose the constants and deformations from the crystal itself, in its own orientation.

    spglib finds the crystal's symmetry from POSCAR (lattice and atoms, as
    given, without rotating anything). The independent elastic constants for
    that orientation and a minimal set of deformations that determine them are
    derived from it, the deformation folders are created, and the setup is saved
    to elastic_vasp_auto.json for auto_post().
    """
    _check_poscar()
    atoms = poscar.read("POSCAR")[0]
    spec, basis = auto.from_structure(atoms, symprec)
    print("Crystal symmetry (spglib, symprec={} A): {}".format(symprec, spec.name.split(": ", 1)[1]))
    print("{} independent elastic constants in this orientation: {}".format(
        len(spec.constants), ", ".join(spec.constants)))
    simple = auto.is_simple(spec)
    if simple:
        print("Stiffness matrix pattern:")
        for row in auto.matrix_pattern(spec):
            print("    " + "  ".join("{:<12}".format(entry) for entry in row))
    else:
        print("In this orientation most entries are linear combinations of these constants; "
              "the full relations are saved in {}.".format(auto.SETUP_FILE))
    print("Deformations ({}):".format(len(spec.deformations)))
    for k, d in enumerate(spec.deformations):
        line = "    {:<8} {}".format(d.folder, d.describe())
        if simple:
            line = "{:<68} a = {}".format(line, spec.energy_expression(k))
        print(line)
    print("Error amplification: at most {:.1f} GPa in C per GPa error in a fitted "
          "coefficient.".format(auto.amplification(spec)))
    generate(spec, overwrite, strains)
    auto.save(auto.SETUP_FILE, spec, basis, symprec, auto.standard_setting(atoms, symprec))
    print("Setup saved to {}; run auto_post() after the VASP runs.".format(auto.SETUP_FILE))
    return spec


def cubic_pre(overwrite=False, strains=None):
    """Cubic: write the strained POSCARs.

    For space groups 195-230 (Laue class m-3m).

    3 independent constants (C11, C12, C44), 3 deformations. Required orientation:
    cube axes along x, y, z. The crystal symmetry is checked with spglib; a warning
    names a working alternative if these deformations cannot measure the crystal as
    oriented.

    Args:
        overwrite: replace existing deformation folders (deletes their contents).
        strains: strain values; default: read strain.dat.
    """
    _check_poscar()
    _check_table("cubic")
    generate("cubic", overwrite, strains)


def orthogonal_pre(overwrite=False, strains=None):
    """Orthorhombic: write the strained POSCARs.

    For space groups 16-74 (Laue class mmm); also works for tetragonal I and higher
    symmetry.

    9 independent constants (C11, C22, C33, C12, C13, C23, C44, C55, C66), 9
    deformations. Required orientation: a along x, b along y, c along z. The crystal
    symmetry is checked with spglib; a warning names a working alternative if these
    deformations cannot measure the crystal as oriented.

    Args:
        overwrite: replace existing deformation folders (deletes their contents).
        strains: strain values; default: read strain.dat.
    """
    _check_poscar()
    _check_table("orthogonal")
    generate("orthogonal", overwrite, strains)


def tetragonal_1_pre(overwrite=False, strains=None):
    """Tetragonal I: write the strained POSCARs.

    For space groups 89-142 (Laue class 4/mmm).

    6 independent constants (C11, C12, C13, C33, C44, C66), 6 deformations. Required
    orientation: c (4-fold axis) along z, a along x or along [110]. The crystal symmetry is checked
    with spglib; a warning names a working alternative if these deformations cannot
    measure the crystal as oriented.

    Args:
        overwrite: replace existing deformation folders (deletes their contents).
        strains: strain values; default: read strain.dat.
    """
    _check_poscar()
    _check_table("tetragonal_1")
    generate("tetragonal_1", overwrite, strains)


def tetragonal_2_pre(overwrite=False, strains=None):
    """Tetragonal II: write the strained POSCARs.

    For space groups 75-88 (Laue class 4/m), including C16.

    7 independent constants (C11, C12, C13, C33, C44, C66, C16), 7 deformations.
    Required orientation: c (4-fold axis) along z; any in-plane orientation (C16 is
    fitted), so it also works for tetragonal I crystals at any angle. The crystal symmetry
    is checked with spglib; a warning names a working alternative if these
    deformations cannot measure the crystal as oriented.

    Args:
        overwrite: replace existing deformation folders (deletes their contents).
        strains: strain values; default: read strain.dat.
    """
    _check_poscar()
    _check_table("tetragonal_2")
    generate("tetragonal_2", overwrite, strains)


def hexagonal_pre(overwrite=False, strains=None):
    """Hexagonal: write the strained POSCARs.

    For space groups 168-194 (Laue class 6/mmm).

    5 independent constants (C11, C12, C13, C33, C44), 5 deformations. Required
    orientation: c (6-fold axis) along z; any in-plane orientation. The crystal
    symmetry is checked with spglib; a warning names a working alternative if these
    deformations cannot measure the crystal as oriented.

    Args:
        overwrite: replace existing deformation folders (deletes their contents).
        strains: strain values; default: read strain.dat.
    """
    _check_poscar()
    if orientation.check_threefold_z(system="hexagonal"):
        _check_table("hexagonal")
    generate("hexagonal", overwrite, strains)


def trigonal_1_pre(overwrite=False, strains=None):
    """Trigonal I: write the strained POSCARs.

    For space groups 149-167 (Laue class -3m); for the '31m' groups 149, 151, 153,
    157, 159, 162 and 163 use trigonal_2.

    6 independent constants (C11, C12, C13, C33, C44, C14), 6 deformations. Required
    orientation: c (3-fold axis) along z and x along a 2-fold axis. The crystal
    symmetry is checked with spglib; a warning names a working alternative if these
    deformations cannot measure the crystal as oriented.

    Args:
        overwrite: replace existing deformation folders (deletes their contents).
        strains: strain values; default: read strain.dat.
    """
    _check_poscar()
    if orientation.check_threefold_z(system="trigonal"):
        _check_table("trigonal_1")
    generate("trigonal_1", overwrite, strains)


def trigonal_2_pre(overwrite=False, strains=None):
    """Trigonal II: write the strained POSCARs.

    For space groups 143-148 (Laue class -3) and the '31m'-type trigonal I groups;
    fits C14 and C15 together.

    7 independent constants (C11, C12, C13, C33, C44, C14, C15), 7 deformations.
    Required orientation: c (3-fold axis) along z; any in-plane orientation. The
    crystal symmetry is checked with spglib; a warning names a working alternative
    if these deformations cannot measure the crystal as oriented.

    Args:
        overwrite: replace existing deformation folders (deletes their contents).
        strains: strain values; default: read strain.dat.
    """
    _check_poscar()
    if orientation.check_threefold_z(system="trigonal"):
        _check_table("trigonal_2")
    generate("trigonal_2", overwrite, strains)


def monoclinic_pre(overwrite=False, strains=None):
    """Monoclinic: write the strained POSCARs.

    For space groups 3-15 (Laue class 2/m), b-unique setting.

    13 independent constants (C11, C22, C33, C12, C13, C23, C44, C55, C66, C15, C25,
    C35, C46), 13 deformations. Required orientation: unique 2-fold axis b along y,
    a and c in the xz plane. The crystal symmetry is checked with spglib; a warning
    names a working alternative if these deformations cannot measure the crystal as
    oriented.

    Args:
        overwrite: replace existing deformation folders (deletes their contents).
        strains: strain values; default: read strain.dat.
    """
    _check_poscar()
    if orientation.check_monoclinic():
        _check_table("monoclinic")
    generate("monoclinic", overwrite, strains)


def triclinic_pre(overwrite=False, strains=None):
    """Triclinic: write the strained POSCARs.

    For space groups 1-2; all 21 constants, so it works for any crystal.

    21 independent constants (C11, C12, C13, C14, C15, C16, C22, C23, C24, C25, C26,
    C33, C34, C35, C36, C44, C45, C46, C55, C56, C66), 21 deformations. Required
    orientation: any. The crystal symmetry is checked with spglib; a warning names a
    working alternative if these deformations cannot measure the crystal as
    oriented.

    Args:
        overwrite: replace existing deformation folders (deletes their contents).
        strains: strain values; default: read strain.dat.
    """
    _check_poscar()
    _check_table("triclinic")
    generate("triclinic", overwrite, strains)
