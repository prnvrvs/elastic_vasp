"""Orientation checks for the strain patterns in elastic_gen.

The strains are applied along the Cartesian x, y, z axes of the POSCAR
lattice, so some crystal systems need a specific alignment:

* monoclinic: unique (2-fold) axis along y
* hexagonal / trigonal / rhombohedral: 3-fold (or 6-fold) axis along z

The checks below test whether the *lattice* has the required rotation
symmetry about that Cartesian axis. This works for any cell choice
(conventional, primitive, centred). If the check fails, a warning explains
where the axis actually is, a rotated copy is written to POSCAR.reoriented
(pure rotation: fractional coordinates and atom order are unchanged, so
POTCAR/MAGMOM ordering stays valid), and pre-processing continues.
"""
import os
import warnings
from itertools import product
from math import gcd, pi

import numpy as np

from elastic_vasp import poscar

LATTICE_TOL = 1e-3
REORIENTED_NAME = "POSCAR.reoriented"


def _rotation(axis, angle):
    n = np.asarray(axis, float)
    n = n / np.linalg.norm(n)
    k = np.array([[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]])
    return np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k


def _is_lattice_symmetry(lattice, rot, tol=LATTICE_TOL):
    """True if rotating the lattice by `rot` maps it onto itself."""
    m = lattice @ rot.T @ np.linalg.inv(lattice)
    return np.allclose(m, np.round(m), atol=tol)


def _find_axes(lattice, angle, max_index=2):
    """Lattice directions [uvw] about which the lattice has `angle` rotation symmetry."""
    found = []
    for uvw in product(range(-max_index, max_index + 1), repeat=3):
        if uvw == (0, 0, 0) or gcd(gcd(abs(uvw[0]), abs(uvw[1])), abs(uvw[2])) != 1:
            continue
        if next(i for i in uvw if i) < 0:
            continue  # report each axis once, with its first nonzero index positive
        direction = np.array(uvw, float) @ lattice
        n = direction / np.linalg.norm(direction)
        if any(abs(abs(n @ other) - 1) < 1e-6 for _, other in found):
            continue
        if _is_lattice_symmetry(lattice, _rotation(n, angle)):
            found.append((uvw, n))
    return found


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


def _check_axis(poscar_path, system, fold, angle, target, target_name):
    atoms, title = poscar.read(poscar_path)
    lattice = np.array(atoms.cell.array)
    target = np.asarray(target, float)
    if _is_lattice_symmetry(lattice, _rotation(target, angle)):
        return True

    msg = [
        "{} pre-processing expects the {} axis along {}, but the lattice in {} has "
        "no {} symmetry about {}.".format(system, fold, target_name, poscar_path, fold, target_name)
    ]
    axes = _find_axes(lattice, angle)
    if not axes:
        msg.append("No {} axis was found in the lattice either; check that the structure "
                   "really is {}.".format(fold, system))
    else:
        uvw, n = axes[0]
        msg.append("The lattice {} axis is along lattice direction [{} {} {}] = Cartesian "
                   "({:.3f}, {:.3f}, {:.3f}).".format(fold, uvw[0], uvw[1], uvw[2], *n))
        if len(axes) > 1:
            msg.append("(The lattice has {} such axes, so it may have higher symmetry than "
                       "{}.)".format(len(axes), system))
        atoms.set_cell(lattice @ _align(n, target).T, scale_atoms=True)  # pure rotation
        out_path = os.path.join(os.path.dirname(poscar_path), REORIENTED_NAME)
        poscar.write(out_path, atoms, title)
        msg.append("A rotated copy with that axis along {} was written to {} "
                   "(fractional coordinates and atom order unchanged). To use it, copy "
                   "it over POSCAR and rerun with overwrite=True.".format(target_name, out_path))
    msg.append("Continuing with the current orientation: the elastic constants will be "
               "wrong if the axis really is misaligned.")
    warnings.warn(" ".join(msg), UserWarning, stacklevel=4)
    return False


def check_monoclinic(poscar_path="POSCAR"):
    """Warn (and write POSCAR.reoriented) unless the unique 2-fold axis is along y."""
    return _check_axis(poscar_path, "monoclinic", "unique (2-fold)", pi, [0, 1, 0], "y")


def check_threefold_z(poscar_path="POSCAR", system="hexagonal/trigonal"):
    """Warn (and write POSCAR.reoriented) unless the 3-fold/6-fold axis is along z."""
    return _check_axis(poscar_path, system, "3-fold/6-fold (c)", 2 * pi / 3, [0, 0, 1], "z")
