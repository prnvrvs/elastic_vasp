"""Orientation checks for the strain patterns in elastic_gen.

The strains are applied along the Cartesian x, y, z axes of the POSCAR
lattice, so some crystal systems need a specific alignment:

* monoclinic: unique (2-fold) axis along y
* hexagonal / trigonal / rhombohedral: 3-fold (or 6-fold) axis along z

The checks use the *crystal* symmetry found by spglib from the lattice and the
atom positions, so a cell whose lattice looks more symmetric than the crystal
(e.g. a monoclinic crystal in a cell with all angles 90 degrees) is still
judged correctly. If the required axis is not where it should be, a warning
says where it is, a rotated copy is written to POSCAR.reoriented (pure
rotation: fractional coordinates and atom order are unchanged, so POTCAR and
MAGMOM ordering stay valid), and pre-processing continues.

Magnetic order is not taken into account: spglib sees only the atoms.
"""
import os
import warnings
from math import pi

import numpy as np

from elastic_vasp import poscar

SYMPREC = 1e-3          # spglib distance tolerance (Angstrom)
AXIS_TOL = 1e-5         # |cos| tolerance for "this axis is along x/y/z" (~0.26 degrees)
REORIENTED_NAME = "POSCAR.reoriented"


def _rotation(axis, angle):
    n = np.asarray(axis, float)
    n = n / np.linalg.norm(n)
    k = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def _align(n, target):
    """Rotation matrix taking unit vector n onto unit vector target."""
    axis = np.cross(n, target)
    cos_angle = float(np.clip(n @ target, -1.0, 1.0))
    if np.linalg.norm(axis) < 1e-12:
        if cos_angle > 0:
            return np.eye(3)
        perp = np.cross(n, [1.0, 0.0, 0.0])
        if np.linalg.norm(perp) < 1e-6:
            perp = np.cross(n, [0.0, 1.0, 0.0])
        return _rotation(perp, pi)
    return _rotation(axis, np.arccos(cos_angle))


def _spglib_cell(atoms):
    return (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)


def space_group(atoms, symprec=SYMPREC):
    import spglib

    return spglib.get_spacegroup(_spglib_cell(atoms), symprec=symprec) or "unknown"


def symmetry_axes(atoms, order, symprec=SYMPREC):
    """Cartesian unit vectors of the crystal's n-fold rotation axes (n = order).

    Improper operations (mirrors, rotoinversions) count through their proper
    part, so e.g. -3 provides a 3-fold axis and 2/m a 2-fold axis. A 6-fold
    axis also appears as a 3-fold axis. Uses the primitive cell, so supercells
    in any orientation keep their full symmetry.
    """
    from elastic_vasp.auto import crystal_rotations

    try:
        rotations = crystal_rotations(atoms, symprec)
    except ValueError:
        return []
    target_cos = np.cos(2 * pi / order)
    axes = []
    for r in rotations:
        if np.linalg.det(r) < 0:
            r = -r
        if abs((np.trace(r) - 1) / 2 - target_cos) > 1e-6:
            continue
        values, vectors = np.linalg.eig(r)
        axis = np.real(vectors[:, np.argmin(np.abs(values - 1))])
        axis /= np.linalg.norm(axis)
        if axis[np.argmax(np.abs(axis))] < 0:
            axis = -axis
        if not any(abs(abs(axis @ other) - 1) < AXIS_TOL for other in axes):
            axes.append(axis)
    return axes


def _check_axis(poscar_path, system, fold, order, target, target_name, symprec):
    atoms, title = poscar.read(poscar_path)
    target = np.asarray(target, float)
    axes = symmetry_axes(atoms, order, symprec)
    if any(abs(abs(axis @ target) - 1) < AXIS_TOL for axis in axes):
        return True

    group = space_group(atoms, symprec)
    msg = ["{} pre-processing expects the crystal's {} axis along {}, but the structure in "
           "{} (space group {}) has none along {}.".format(
               system, fold, target_name, poscar_path, group, target_name)]
    if not axes:
        msg.append("It has no {} axis at all (spglib, symprec={} A), so it does not look {}; "
                   "check the structure or use a lower-symmetry system such as "
                   "triclinic_pre().".format(fold, symprec, system))
    else:
        n = axes[0]
        msg.append("Its {} axis is along Cartesian ({:.3f}, {:.3f}, {:.3f}).".format(fold, *n))
        if len(axes) > 1:
            msg.append("(It has {} such axes, so the crystal has higher symmetry than {}; "
                       "check that {} is the right system.)".format(len(axes), system, system))
        atoms.set_cell(atoms.cell.array @ _align(n, target).T, scale_atoms=True)  # pure rotation
        out_path = os.path.join(os.path.dirname(poscar_path), REORIENTED_NAME)
        poscar.write(out_path, atoms, title)
        msg.append("A rotated copy with that axis along {} was written to {} "
                   "(fractional coordinates and atom order unchanged). To use it, copy "
                   "it over POSCAR and rerun with overwrite=True.".format(target_name, out_path))
    msg.append("Continuing with the current orientation: the elastic constants will be "
               "wrong if the axis really is misaligned.")
    warnings.warn(" ".join(msg), UserWarning, stacklevel=4)
    return False


def check_monoclinic(poscar_path="POSCAR", symprec=SYMPREC):
    """Warn (and write POSCAR.reoriented) unless the crystal's unique 2-fold axis is along y."""
    return _check_axis(poscar_path, "monoclinic", "2-fold (unique)", 2, [0, 1, 0], "y", symprec)


def check_threefold_z(poscar_path="POSCAR", system="hexagonal/trigonal", symprec=SYMPREC):
    """Warn (and write POSCAR.reoriented) unless the crystal's 3-fold/6-fold axis is along z."""
    return _check_axis(poscar_path, system, "3-fold/6-fold (c)", 3, [0, 0, 1], "z", symprec)


def auto_noise_tol():
    from elastic_vasp.auto import NOISE_TOL
    return NOISE_TOL


def table_covers(spec, basis):
    """True if the named table's layout can represent every symmetry-allowed stiffness matrix."""
    n = len(spec.constants)
    table = np.array([spec.matrix(np.eye(n)[j]).ravel() for j in range(n)]).T
    for b in basis.values():
        target = np.asarray(b).ravel()
        coef, *_ = np.linalg.lstsq(table, target, rcond=None)
        if np.abs(table @ coef - target).max() > auto_noise_tol():
            return False
    return True


def check_table(system, poscar_path="POSCAR", symprec=SYMPREC):
    """Warn unless the named system's deformations can measure this crystal's stiffness matrix.

    The symmetry-allowed stiffness matrices of the crystal (from spglib, in the
    POSCAR's own orientation) must all be representable by the table's layout;
    otherwise some constants are never measured (e.g. C16 of a tetragonal II
    crystal with orthogonal_pre) or are tied together wrongly.
    """
    from elastic_vasp import auto
    from elastic_vasp.systems import SYSTEMS

    atoms = poscar.read(poscar_path)[0]
    basis = auto.invariant_basis(auto.crystal_rotations(atoms, symprec))
    spec = SYSTEMS[system]
    if table_covers(spec, basis):
        return True

    n = len(spec.constants)
    in_table = np.zeros((6, 6), bool)
    for j in range(n):
        in_table |= np.abs(spec.matrix(np.eye(n)[j])) > 1e-9
    allowed = np.zeros((6, 6), bool)
    for b in basis.values():
        allowed |= np.abs(np.asarray(b)) > auto_noise_tol()
    missing = ["C{}{}".format(i + 1, j + 1) for i in range(6) for j in range(i, 6)
               if allowed[i, j] and not in_table[i, j]]
    working = sorted((name for name, table in SYSTEMS.items() if table_covers(table, basis)),
                     key=lambda name: (len(SYSTEMS[name].constants), name))
    msg = ["{}_pre() cannot measure this crystal's stiffness matrix: {} (space group {}) has {} "
           "independent elastic constants in this orientation".format(
               system, poscar_path, space_group(atoms, symprec), len(basis))]
    if missing:
        msg[0] += ", and {} can be non-zero, which the {} deformations assume to be zero.".format(
            ", ".join(missing), system)
    else:
        msg[0] += ", related differently than the {} deformations assume (e.g. the crystal axes are " \
                  "not where they expect).".format(system)
    alternatives = ["{}_pre()".format(name) for name in working[:2]]
    msg.append("Use {}auto_pre(), which works in any orientation.".format(
        " or ".join(alternatives) + " or " if alternatives else ""))
    msg.append("Continuing: the elastic constants will be wrong.")
    warnings.warn(" ".join(msg), UserWarning, stacklevel=4)
    return False
