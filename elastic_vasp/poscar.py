"""Minimal POSCAR helpers shared by pre-processing, volume reading and the
orientation checks.

VASP's second POSCAR line is the scaling factor:

* one positive value s: lattice vectors are multiplied by s
* one negative value -V: lattice is rescaled so the cell volume is V
* three values (VASP 6): x, y and z components are scaled separately
"""
import numpy as np


def scale_values(line):
    """Numbers on the scale line (ignores trailing comments)."""
    values = []
    for token in line.split():
        try:
            values.append(float(token))
        except ValueError:
            break
        if len(values) == 3:
            break
    if not values:
        raise ValueError("POSCAR scale line has no number: {!r}".format(line))
    if len(values) == 2:
        raise ValueError("POSCAR scale line must have 1 or 3 values: {!r}".format(line))
    return values


def raw_lattice(lines):
    """The three lattice-vector rows exactly as written (before scaling)."""
    if len(lines) < 5:
        raise ValueError("POSCAR is too short")
    return np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)])


def is_plain_scale(lines):
    """True if the scale line is a single positive number."""
    values = scale_values(lines[1])
    return len(values) == 1 and values[0] > 0


def lattice(lines):
    """Actual Cartesian lattice (rows = vectors, Angstrom) with the scale applied."""
    raw = raw_lattice(lines)
    values = scale_values(lines[1])
    if len(values) == 3:
        if min(values) <= 0:
            raise ValueError("POSCAR per-axis scale factors must be positive: {}".format(values))
        return raw * np.array(values)
    scale = values[0]
    if scale == 0:
        raise ValueError("POSCAR scale factor is zero")
    if scale < 0:
        return raw * (abs(scale) / abs(np.linalg.det(raw))) ** (1.0 / 3.0)
    return raw * scale


def read(path):
    with open(path) as fh:
        return fh.readlines()
