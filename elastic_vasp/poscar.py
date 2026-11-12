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

    with open(path) as fh:
        lines = fh.readlines()
    title = lines[0].rstrip("\n") if lines else ""
    source = str(path)
    species = _species_tokens(lines)
    if species is not None:
        elements = [_element(token) for token in species]
        if elements != species:
            # POTCAR-style labels (Fe_pv, O_s, Fe/<hash>) are fine for VASP, which
            # takes species from POTCAR, but ASE needs plain element symbols.
            lines = lines[:5] + ["  " + "  ".join(elements) + "\n"] + lines[6:]
            source = io.StringIO("".join(lines))
    try:
        atoms = read_vasp(source)
    except Exception as exc:  # ASE raises many types for malformed files
        raise ValueError("{} could not be read by ASE ({}: {})".format(
            path, type(exc).__name__, exc)) from None
    return atoms, title


def _species_tokens(lines):
    """Tokens of the VASP 5 species line, or None for VASP 4 files (no such line)."""
    if len(lines) < 6:
        return None
    tokens = lines[5].split()
    if not tokens:
        return None
    try:
        int(tokens[0])
        return None  # line 6 already holds the atom counts
    except ValueError:
        return tokens


def _element(label):
    """'Fe_pv' -> 'Fe', 'O_s' -> 'O', 'Fe/abc123' -> 'Fe'."""
    return label.split("/", 1)[0].split("_", 1)[0]


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
