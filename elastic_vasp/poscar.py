"""POSCAR reading and writing, done by ASE.

ASE handles every POSCAR form VASP accepts that matters here: a positive,
negative (target volume) or per-axis (VASP 6) scale factor, Direct or
Cartesian coordinates, and selective dynamics. Structures are always written
with scale 1.0 and Direct coordinates, keeping the original title line and the
original atom order.
"""
import io

import numpy as np

# Older ASE silently ignores negative (target-volume) and per-axis scale factors.
MIN_ASE = (3, 23)


def _check_ase():
    import ase

    version = tuple(int(x) for x in ase.__version__.split(".")[:2] if x.isdigit())
    if version < MIN_ASE:
        raise ImportError(
            "elastic_vasp needs ASE >= {}.{} (found {}); older versions mis-read POSCAR "
            "scale factors. Upgrade with: pip install -U ase".format(*MIN_ASE, ase.__version__))


def read(path):
    """(ase.Atoms, title line) from a POSCAR file.

    Raises ValueError with a readable message if ASE cannot parse the file.
    """
    _check_ase()
    from ase.io.vasp import read_vasp

    try:
        atoms = read_vasp(str(path))
    except Exception as exc:  # ASE raises many types for malformed files
        raise ValueError("{} could not be read by ASE ({}: {})".format(
            path, type(exc).__name__, exc)) from None
    with open(path) as fh:
        title = fh.readline().rstrip("\n")
    return atoms, title


def write(path, atoms, title="POSCAR"):
    """Write atoms as a VASP 5 POSCAR: scale 1.0, Direct coordinates, original order."""
    from ase.io.vasp import write_vasp

    buffer = io.StringIO()
    write_vasp(buffer, atoms, direct=True, sort=False)
    lines = buffer.getvalue().splitlines(keepends=True)
    lines[0] = title + "\n"
    with open(path, "w") as fh:
        fh.writelines(lines)


def lattice(path):
    """Actual Cartesian lattice (rows = vectors, Angstrom), scale factor applied."""
    return np.array(read(path)[0].cell.array)


def volume(path):
    return abs(float(np.linalg.det(lattice(path))))
