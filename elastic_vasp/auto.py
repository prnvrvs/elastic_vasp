"""Orientation-aware setup: constants and deformations from the crystal itself.

Instead of assuming a standard orientation (c along z, unique axis along y,
...), the crystal's symmetry operations are taken from spglib in the user's own
Cartesian frame. From them:

1. The stiffness matrices invariant under every operation form a linear
   subspace. Its reduced row-echelon basis over the entries
   C11, C12, ..., C16, C22, ..., C66 (row-major upper triangle) names the
   independent constants: for standard orientations these are exactly the
   textbook sets (cubic C11, C12, C44; hexagonal C11, C12, C13, C33, C44; ...),
   for other orientations they are the entries that are independent there.
2. A minimal set of deformations is picked greedily from simple candidates
   (one strain component, or two with equal or opposite sign) so that the
   energies determine every independent constant, keeping the linear system as
   well conditioned as possible.

The result is an ordinary systems.CrystalSystem, so pre- and post-processing
work exactly as for the named systems. It is saved to JSON by pre-processing
so post-processing uses exactly the deformations that were generated.
"""
import json
from itertools import combinations

import numpy as np

from elastic_vasp.systems import CrystalSystem, Deformation

SETUP_FILE = "elastic_vasp_auto.json"
SYMPREC = 1e-3
_TOL = 1e-8
# Relaxed structures are symmetric only to ~1e-7, which leaves entries of that
# size where the ideal crystal has zeros. Comparisons and displays treat
# anything below NOISE_TOL (relative to pivots of 1) as zero.
NOISE_TOL = 1e-6

# Voigt index -> Cartesian index pair
_VOIGT = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]
# independent entries of a symmetric 6x6 matrix, row-major upper triangle
_UPPER = [(i, j) for i in range(6) for j in range(i, 6)]
_NAMES = ["C{}{}".format(i + 1, j + 1) for i, j in _UPPER]


def bond_matrix(rotation):
    """6x6 matrix M with C' = M C M^T for a Cartesian rotation (engineering shear strains)."""
    r = np.asarray(rotation, float)
    m = np.zeros((6, 6))
    for a, (i, j) in enumerate(_VOIGT):
        for b, (k, l) in enumerate(_VOIGT):
            m[a, b] = r[i, k] * r[j, l] if b < 3 else r[i, k] * r[j, l] + r[i, l] * r[j, k]
    return m


def _to_matrix(vector):
    c = np.zeros((6, 6))
    for value, (i, j) in zip(vector, _UPPER):
        c[i, j] = c[j, i] = value
    return c


def _to_vector(matrix):
    return np.array([matrix[i, j] for i, j in _UPPER])


def primitive_cell(atoms, symprec=SYMPREC):
    """spglib cell (lattice, positions, numbers) of the primitive cell, same orientation.

    For a supercell (e.g. an fcc crystal in a cell with x || [1-10], y || [11-2],
    z || [111]) spglib only returns the operations that also map that particular
    supercell onto itself, which can be a small subgroup of the crystal's point
    group. The primitive cell, without idealization (no rotation), has them all.
    """
    import spglib

    cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
    primitive = spglib.standardize_cell(cell, to_primitive=True, no_idealize=True, symprec=symprec)
    return primitive if primitive is not None else cell


def crystal_rotations(atoms, symprec=SYMPREC):
    """Cartesian rotation matrices of the crystal's point group (from spglib)."""
    import spglib

    cell = primitive_cell(atoms, symprec)
    data = spglib.get_symmetry(cell, symprec=symprec)
    if not data:
        raise ValueError("spglib could not determine the symmetry (symprec={})".format(symprec))
    lattice = np.asarray(cell[0])
    to_cart, from_cart = lattice.T, np.linalg.inv(lattice.T)
    rotations = []
    for w in data["rotations"]:
        r = to_cart @ w @ from_cart
        # A relaxed cell is only symmetric to ~1e-7, so r is not exactly orthogonal;
        # use the nearest exact rotation (or rotoinversion) instead.
        u, _, vt = np.linalg.svd(r)
        r = u @ vt
        if not any(np.allclose(r, other, atol=1e-6) for other in rotations):
            rotations.append(r)
    return rotations


def standard_setting(atoms, symprec=SYMPREC):
    """(rotation, space-group number) taking this crystal to its standard setting.

    The rotation R (from spglib's ``std_rotation_matrix``) maps Cartesian vectors
    of the input frame to the frame of the crystal's standard (idealized)
    setting: cube axes along x, y, z; hexagonal / trigonal c along z with a
    along x; monoclinic unique axis b along y; and so on. A stiffness matrix
    then transforms as ``to_standard(C, R)``.
    """
    import spglib

    dataset = spglib.get_symmetry_dataset(primitive_cell(atoms, symprec), symprec=symprec)
    if dataset is None:
        raise ValueError("spglib could not determine the symmetry (symprec={})".format(symprec))
    get = (lambda key: getattr(dataset, key)) if hasattr(dataset, "number") else dataset.__getitem__
    rotation = np.array(get("std_rotation_matrix"), float)
    u, _, vt = np.linalg.svd(rotation)
    return u @ vt, int(get("number"))


def to_standard(c, rotation):
    """Stiffness matrix C (input frame) expressed in the standard frame given by `rotation`."""
    m = bond_matrix(rotation)
    return m @ np.asarray(c, float) @ m.T


def standard_constants(c_standard, number):
    """Named workflow and {constant: value} for a stiffness matrix in the standard setting.

    Returns (name, values, residual): the textbook constants of the crystal
    system that the space group belongs to, and how far (GPa) the matrix is
    from that system's form (should be ~0; large values mean the symmetry and
    the matrix disagree).
    """
    from elastic_vasp.symmetry import classify
    from elastic_vasp.systems import SYSTEMS

    name = classify(number)[2]
    table = SYSTEMS[name]
    n = len(table.constants)
    design = np.array([table.matrix(np.eye(n)[j]).ravel() for j in range(n)]).T
    values, *_ = np.linalg.lstsq(design, np.asarray(c_standard).ravel(), rcond=None)
    residual = float(np.abs(design @ values - np.asarray(c_standard).ravel()).max())
    return name, dict(zip(table.constants, values)), residual


def invariant_basis(rotations):
    """{name: 6x6 matrix}: reduced row-echelon basis of the symmetry-allowed stiffness matrices."""
    eye = np.eye(len(_UPPER))
    blocks = []
    for r in rotations:
        m = bond_matrix(r)
        t = np.array([_to_vector(m @ _to_matrix(e) @ m.T) for e in eye]).T  # 21x21
        blocks.append(t - eye)
    a = np.vstack(blocks) if blocks else np.zeros((1, len(_UPPER)))
    _, s, vt = np.linalg.svd(a)
    s = np.concatenate([s, np.zeros(vt.shape[0] - len(s))])
    null = vt[s < 1e-6 * max(1.0, s.max())]
    rows = _rref(null)
    basis = {}
    for row in rows:
        pivot = int(np.argmax(np.abs(row) > _TOL))
        row = np.where(np.abs(row) < 1e-10, 0.0, row)
        basis[_NAMES[pivot]] = _to_matrix(row)
    return basis


def _rref(m):
    m = np.array(m, float)
    rows, cols = m.shape
    lead = 0
    for c in range(cols):
        if lead >= rows:
            break
        p = lead + int(np.argmax(np.abs(m[lead:, c])))
        if abs(m[p, c]) < _TOL:
            continue
        m[[lead, p]] = m[[p, lead]]
        m[lead] /= m[lead, c]
        for r in range(rows):
            if r != lead:
                m[r] -= m[r, c] * m[lead]
        lead += 1
    return m[:lead]


def _candidates():
    """Simple deformations in order of preference, with volume-conserving corrections."""
    out = []
    for i in range(6):
        e = [0] * 6
        e[i] = 1
        correction = (i - 3, 4) if i >= 3 else None  # shear: keep volume via e(i-3)
        out.append((e, correction))
    for i, j in combinations(range(6), 2):
        for sign in (1, -1):
            e = [0] * 6
            e[i], e[j] = 1, sign
            correction = None
            if sign < 0 and j < 3:  # (s, -s) on two normal strains: fix the third
                correction = (3 - i - j, 1)
            out.append((e, correction))
    return out


def _folder_name(strain):
    parts = []
    for index, value in enumerate(strain):
        if value:
            label = "e{}".format(index + 1)
            parts.append(label if not parts else ("+" if value > 0 else "-") + label)
    return "".join(parts)


def _energy_row(strain, basis_matrices):
    e = np.asarray(strain, float)
    return np.array([0.5 * e @ b @ e for b in basis_matrices])


def _orthonormal(basis):
    """Orthonormal (Frobenius) basis of the same subspace of stiffness matrices."""
    vectors = np.array([np.asarray(m).ravel() for m in basis.values()]).T
    q, _ = np.linalg.qr(vectors)
    return [q[:, j].reshape(6, 6) for j in range(q.shape[1])]


def choose_deformations(basis):
    """Minimal set of deformations that determines the stiffness matrix as robustly as possible.

    With an orthonormal basis of the allowed matrices, the error in C caused by
    errors in the fitted energy coefficients is at most 1/sigma_min of the
    linear system, independent of how the constants are named. Deformations
    are added greedily to keep sigma_min as large as possible (rows divided
    by e.e so one- and two-component strains compare fairly).
    """
    matrices = _orthonormal(basis)
    n = len(matrices)
    chosen, rows = [], []
    pool = _candidates()
    while len(rows) < n:
        best = None
        for strain, correction in pool:
            if any(strain == c[0] for c in chosen):
                continue
            e = np.asarray(strain, float)
            row = _energy_row(strain, matrices) / (e @ e)
            if np.linalg.norm(row) < _TOL:
                continue
            s = np.linalg.svd(np.vstack(rows + [row]), compute_uv=False)
            rank = int(np.sum(s > 1e-8))
            score = (rank, s[rank - 1] if rank else 0.0)
            if best is None or score[0] > best[0][0] or (
                    score[0] == best[0][0] and score[1] > best[0][1] * (1 + 1e-9)):
                best = (score, strain, correction, row)
        if best is None or best[0][0] <= len(rows):
            raise ValueError("could not find deformations that determine all constants")
        chosen.append((best[1], best[2]))
        rows.append(best[3])
    return tuple(Deformation(_folder_name(s), tuple(s), c) for s, c in chosen)


def amplification(spec):
    """Worst-case GPa error in C per GPa error in the fitted coefficients (spectral norm)."""
    design_inverse = np.linalg.pinv(spec.design_matrix())
    n = len(spec.constants)
    basis = np.array([spec.matrix(np.eye(n)[j]).ravel() for j in range(n)]).T
    return float(np.linalg.norm(basis @ design_inverse, 2))


def _layout(basis):
    names = list(basis)

    def layout(ns):
        return sum(getattr(ns, name) * basis[name] for name in names).tolist()
    return layout


def from_structure(atoms, symprec=SYMPREC):
    """CrystalSystem for this crystal in its own orientation."""
    from elastic_vasp.orientation import space_group

    basis = invariant_basis(crystal_rotations(atoms, symprec))
    spec = CrystalSystem(
        "auto: {}".format(space_group(atoms, symprec)), tuple(basis), _layout(basis),
        choose_deformations(basis))
    return spec, basis


def matrix_pattern(spec):
    """6x6 grid of strings showing each entry in terms of the independent constants."""
    grid = []
    for i in range(6):
        row = []
        for j in range(6):
            terms = []
            for name in spec.constants:
                coeff = spec.matrix([1.0 if n == name else 0.0 for n in spec.constants])[i, j]
                if abs(coeff) > NOISE_TOL:
                    terms.append((coeff, name))
            row.append(_format_terms(terms))
        grid.append(row)
    return grid


def is_simple(spec):
    """True if every entry is 0 or a combination with coefficients +-1, +-0.5 or +-2."""
    n = len(spec.constants)
    for j in range(n):
        m = spec.matrix(np.eye(n)[j])
        values = np.abs(m[np.abs(m) > NOISE_TOL])
        if not all(min(abs(v - t) for t in (0.5, 1.0, 2.0)) < NOISE_TOL for v in values):
            return False
    return True


def _format_terms(terms):
    if not terms:
        return "0"
    out = ""
    for k, (coeff, name) in enumerate(terms):
        sign = "-" if coeff < 0 else "+"
        mag = abs(coeff)
        body = name if abs(mag - 1) < 1e-9 else "{:.4g}*{}".format(mag, name)
        out += ("-" if sign == "-" else "") + body if k == 0 else " {} {}".format(sign, body)
    return out


def save(path, spec, basis, symprec, standard=None):
    """Write the setup; `standard` = (rotation, space-group number) from standard_setting()."""
    data = {
        "format": 1,
        "system": spec.name,
        "symprec": symprec,
        "standard_rotation": np.asarray(standard[0]).tolist() if standard else None,
        "space_group_number": standard[1] if standard else None,
        "constants": list(spec.constants),
        "basis": {name: np.asarray(m).tolist() for name, m in basis.items()},
        "deformations": [
            {"folder": d.folder, "strain": list(d.strain),
             "correction": list(d.correction) if d.correction else None}
            for d in spec.deformations],
    }
    with open(path, "w") as fh:
        json.dump(data, fh, indent=1)


def load_standard(path):
    """(rotation, space-group number) saved by auto_pre, or None for older setup files."""
    with open(path) as fh:
        data = json.load(fh)
    if not data.get("standard_rotation"):
        return None
    return np.array(data["standard_rotation"], float), int(data["space_group_number"])


def load(path):
    with open(path) as fh:
        data = json.load(fh)
    basis = {name: np.array(m, float) for name, m in data["basis"].items()}
    deformations = tuple(
        Deformation(d["folder"], tuple(d["strain"]),
                    tuple(d["correction"]) if d["correction"] else None)
        for d in data["deformations"])
    return CrystalSystem(data["system"], tuple(data["constants"]), _layout(basis), deformations)
