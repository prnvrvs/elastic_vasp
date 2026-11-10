"""Crystal-system tables: one definition drives both pre- and post-processing.

Each crystal system lists

* its independent elastic constants,
* how they fill the 6x6 Voigt stiffness matrix, and
* the deformations: folder name, linear strain direction e (Voigt, with
  engineering shear strains) and an optional second-order correction term.

For a deformation of size s, the energy density is E/V = a s^2 with
a = 1/2 e^T C e (to leading order). ``coeff()`` fits a for every folder, and
because C is linear in the independent constants, ``solve()`` recovers them
from a small linear system. Pre-processing writes exactly the strains listed
here, so the two sides can no longer disagree.
"""
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np


def _square_ratio(value, denominator):
    if abs(denominator) < 1e-15:
        raise ValueError("Invalid strain value {} produces division by zero".format(value))
    return (value ** 2) / denominator


@dataclass(frozen=True)
class Deformation:
    folder: str
    strain: tuple
    # (component index, k): adds s^2 / (k - s^2) to that component, a
    # second-order term that keeps the volume (nearly) constant
    correction: tuple = None

    def at(self, s):
        """Voigt strain (e1..e6) applied for deformation size s."""
        values = [0.0 if c == 0 else c * s for c in self.strain]
        if self.correction is not None:
            index, k = self.correction
            values[index] = _square_ratio(s, k - s ** 2)
        return tuple(values)

    def describe(self):
        text = "(" + ", ".join("%g" % c for c in self.strain) + ") * s"
        if self.correction is not None:
            index, k = self.correction
            text += ", e{} = s^2/({:g}-s^2)".format(index + 1, k)
        return text


@dataclass(frozen=True)
class CrystalSystem:
    name: str
    constants: tuple
    layout: object  # callable: namespace of constants -> 6x6 nested list
    deformations: tuple
    stability_text: str = "stiffness matrix is positive definite (all eigenvalues positive)"

    def matrix(self, values):
        """6x6 stiffness matrix from a sequence of independent constants."""
        ns = SimpleNamespace(**dict(zip(self.constants, values)))
        return np.array(self.layout(ns), dtype=float)

    def design_matrix(self):
        """A with a = A p, where p are the independent constants."""
        basis = np.eye(len(self.constants))
        unit_matrices = [self.matrix(row) for row in basis]
        design = np.empty((len(self.deformations), len(self.constants)))
        for k, deformation in enumerate(self.deformations):
            e = np.array(deformation.strain, dtype=float)
            for j, cj in enumerate(unit_matrices):
                design[k, j] = 0.5 * e @ cj @ e
        return design

    def energy_expression(self, k, design=None):
        """Human-readable a_k = 1/2 e^T C e for deformation k.

        Pass a precomputed design_matrix() when calling this for every row.
        """
        if design is None:
            design = self.design_matrix()
        terms = []
        for coeff, name in zip(design[k], self.constants):
            if abs(coeff) < 1e-12:
                continue
            sign = "-" if coeff < 0 else "+"
            mag = abs(coeff)
            body = name if abs(mag - 1) < 1e-12 else "{:g}*{}".format(mag, name)
            terms.append((sign, body))
        if not terms:
            return "0"
        first_sign, first = terms[0]
        text = ("-" if first_sign == "-" else "") + first
        for sign, body in terms[1:]:
            text += " {} {}".format(sign, body)
        return text

    def solve(self, coefficients):
        """Independent constants and the 6x6 matrix from fitted coefficients.

        `coefficients` maps folder name -> fitted a (GPa).
        """
        design = self.design_matrix()
        rank = np.linalg.matrix_rank(design)
        if rank < len(self.constants):
            raise ValueError(
                "{}: the deformations determine only {} of {} constants".format(
                    self.name, rank, len(self.constants)))
        rhs = np.array([coefficients[d.folder] for d in self.deformations], dtype=float)
        values, *_ = np.linalg.lstsq(design, rhs, rcond=None)
        return dict(zip(self.constants, values)), self.matrix(values)


def _d(folder, strain, correction=None):
    return Deformation(folder, tuple(strain), correction)


# --- stiffness layouts --------------------------------------------------------

def _cubic(c):
    return [[c.C11, c.C12, c.C12, 0, 0, 0],
            [c.C12, c.C11, c.C12, 0, 0, 0],
            [c.C12, c.C12, c.C11, 0, 0, 0],
            [0, 0, 0, c.C44, 0, 0],
            [0, 0, 0, 0, c.C44, 0],
            [0, 0, 0, 0, 0, c.C44]]


def _hexagonal(c):
    c66 = (c.C11 - c.C12) / 2
    return [[c.C11, c.C12, c.C13, 0, 0, 0],
            [c.C12, c.C11, c.C13, 0, 0, 0],
            [c.C13, c.C13, c.C33, 0, 0, 0],
            [0, 0, 0, c.C44, 0, 0],
            [0, 0, 0, 0, c.C44, 0],
            [0, 0, 0, 0, 0, c66]]


def _orthorhombic(c):
    return [[c.C11, c.C12, c.C13, 0, 0, 0],
            [c.C12, c.C22, c.C23, 0, 0, 0],
            [c.C13, c.C23, c.C33, 0, 0, 0],
            [0, 0, 0, c.C44, 0, 0],
            [0, 0, 0, 0, c.C55, 0],
            [0, 0, 0, 0, 0, c.C66]]


def _trigonal_1(c):
    # point groups 32, 3m, -3m with x along a 2-fold axis
    c66 = (c.C11 - c.C12) / 2
    return [[c.C11, c.C12, c.C13, c.C14, 0, 0],
            [c.C12, c.C11, c.C13, -c.C14, 0, 0],
            [c.C13, c.C13, c.C33, 0, 0, 0],
            [c.C14, -c.C14, 0, c.C44, 0, 0],
            [0, 0, 0, 0, c.C44, c.C14],
            [0, 0, 0, 0, c.C14, c66]]


def _trigonal_2(c):
    # point groups 3, -3: adds C15 with C25 = -C15, C46 = -C15
    c66 = (c.C11 - c.C12) / 2
    return [[c.C11, c.C12, c.C13, c.C14, c.C15, 0],
            [c.C12, c.C11, c.C13, -c.C14, -c.C15, 0],
            [c.C13, c.C13, c.C33, 0, 0, 0],
            [c.C14, -c.C14, 0, c.C44, 0, -c.C15],
            [c.C15, -c.C15, 0, 0, c.C44, c.C14],
            [0, 0, 0, -c.C15, c.C14, c66]]


def _monoclinic(c):
    # unique axis b along y
    return [[c.C11, c.C12, c.C13, 0, c.C15, 0],
            [c.C12, c.C22, c.C23, 0, c.C25, 0],
            [c.C13, c.C23, c.C33, 0, c.C35, 0],
            [0, 0, 0, c.C44, 0, c.C46],
            [c.C15, c.C25, c.C35, 0, c.C55, 0],
            [0, 0, 0, c.C46, 0, c.C66]]


_TRICLINIC_NAMES = tuple("C{}{}".format(i, j) for i in range(1, 7) for j in range(i, 7))


def _triclinic(c):
    values = vars(c)
    return [[values["C{}{}".format(min(i, j), max(i, j))] for j in range(1, 7)]
            for i in range(1, 7)]


# --- deformation sets ---------------------------------------------------------
# Folder names, strain directions and corrections are identical to the
# original hand-written versions, so existing calculation folders still work.

_ORTHORHOMBIC_NORMAL = (
    _d("C11_C22_C12_I", (1, -1, 0, 0, 0, 0), (2, 1)),
    _d("C11_C22_C12_II", (1, 1, 0, 0, 0, 0)),
    _d("C22_C33_C23_I", (0, 1, -1, 0, 0, 0), (0, 1)),
    _d("C22_C33_C23_II", (0, 1, 1, 0, 0, 0)),
    _d("C11_C33_C13_I", (1, 0, -1, 0, 0, 0), (1, 1)),
    _d("C11_C33_C13_II", (1, 0, 1, 0, 0, 0)),
    _d("C44", (0, 0, 0, 1, 0, 0), (0, 4)),
)

_TRIGONAL_COMMON = (
    _d("C11_C12_I", (1, 1, 0, 0, 0, 0)),
    _d("C11_C12_II", (1, -1, 0, 0, 0, 0), (2, 1)),
    _d("C11_C33_C13_I", (1, 0, 1, 0, 0, 0)),
    _d("C11_C33_C13_II", (1, 0, -1, 0, 0, 0), (1, 1)),
    _d("C11_C44_C14_I", (1, 0, 0, 1, 0, 0)),
    _d("C11_C44_C14_II", (1, 0, 0, -1, 0, 0)),
)

SYSTEMS = {
    "cubic": CrystalSystem(
        "cubic", ("C11", "C12", "C44"), _cubic,
        (
            _d("C11_C12_I", (1, -1, 0, 0, 0, 0)),
            _d("C11_C12_II", (1, 1, 1, 0, 0, 0)),
            _d("C44", (0, 0, 0, 1, 1, 1)),
        ),
        "C11-C12 > 0 ; C11+2C12 > 0 ; C44 > 0 (equivalent to positive definiteness)",
    ),
    "orthogonal": CrystalSystem(
        "orthogonal", ("C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66"),
        _orthorhombic,
        _ORTHORHOMBIC_NORMAL + (
            _d("C55", (0, 0, 0, 0, 1, 0), (1, 4)),
            _d("C66", (0, 0, 0, 0, 0, 1), (2, 4)),
        ),
    ),
    "hexagonal": CrystalSystem(
        "hexagonal", ("C11", "C12", "C13", "C33", "C44"), _hexagonal,
        (
            _d("C11_C12_I", (1, 1, 0, 0, 0, 0)),
            _d("C11_C12_II", (1, -1, 0, 0, 0, 0), (2, 1)),
            _d("C11_C33_C13_I", (1, 0, -1, 0, 0, 0), (1, 1)),
            _d("C44", (0, 0, 0, 1, 1, 0)),
            _d("C11_C12_C13_C33", (1, 1, 1, 0, 0, 0)),
        ),
        "C11 > |C12| ; 2C13^2 < C33(C11+C12) ; C44 > 0 ; C66 > 0 "
        "(equivalent to positive definiteness)",
    ),
    "trigonal_1": CrystalSystem(
        "trigonal_1", ("C11", "C12", "C13", "C33", "C44", "C14"), _trigonal_1,
        _TRIGONAL_COMMON,
    ),
    "trigonal_2": CrystalSystem(
        "trigonal_2", ("C11", "C12", "C13", "C33", "C44", "C14", "C15"), _trigonal_2,
        _TRIGONAL_COMMON + (_d("C11_C44_C15", (1, 0, 0, 0, -1, 0)),),
    ),
    "monoclinic": CrystalSystem(
        "monoclinic",
        ("C11", "C22", "C33", "C12", "C13", "C23", "C44", "C55", "C66",
         "C15", "C25", "C35", "C46"),
        _monoclinic,
        _ORTHORHOMBIC_NORMAL + (
            _d("C11_C55_C15_I", (1, 0, 0, 0, 1, 0)),
            _d("C11_C55_C15_II", (1, 0, 0, 0, -1, 0)),
            _d("C22_C55_C25", (0, 1, 0, 0, 1, 0)),
            _d("C33_C55_C35", (0, 0, 1, 0, 1, 0)),
            _d("C44_C66_C46_I", (0, 0, 0, 1, 0, 1)),
            _d("C44_C66_C46_II", (0, 0, 0, 1, 0, -1)),
        ),
    ),
    "triclinic": CrystalSystem(
        "triclinic", _TRICLINIC_NAMES, _triclinic,
        (
            _d("C11", (1, 0, 0, 0, 0, 0)),
            _d("C11_C22_C12", (1, 1, 0, 0, 0, 0)),
            _d("C22", (0, 1, 0, 0, 0, 0)),
            _d("C22_C33_C23", (0, 1, 1, 0, 0, 0)),
            _d("C33", (0, 0, 1, 0, 0, 0)),
            _d("C11_C33_C13", (1, 0, 1, 0, 0, 0)),
            _d("C44", (0, 0, 0, 1, 0, 0)),
            _d("C55", (0, 0, 0, 0, 1, 0)),
            _d("C66", (0, 0, 0, 0, 0, 1)),
            _d("C11_C44_C14", (1, 0, 0, 1, 0, 0)),
            _d("C11_C55_C15", (1, 0, 0, 0, 1, 0)),
            _d("C11_C66_C16", (1, 0, 0, 0, 0, 1)),
            _d("C22_C44_C24", (0, 1, 0, 1, 0, 0)),
            _d("C22_C55_C25", (0, 1, 0, 0, 1, 0)),
            _d("C22_C66_C26", (0, 1, 0, 0, 0, 1)),
            _d("C33_C44_C34", (0, 0, 1, 1, 0, 0)),
            _d("C33_C55_C35", (0, 0, 1, 0, 1, 0)),
            _d("C33_C66_C36", (0, 0, 1, 0, 0, 1)),
            _d("C44_C55_C45", (0, 0, 0, 1, 1, 0)),
            _d("C44_C66_C46", (0, 0, 0, 1, 0, 1)),
            _d("C55_C66_C56", (0, 0, 0, 0, 1, 1)),
        ),
    ),
}
