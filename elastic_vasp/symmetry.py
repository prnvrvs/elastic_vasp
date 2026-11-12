"""Symmetry report for a structure: space group, constants and which workflow to use.

``check_symmetry()`` only reads the POSCAR; it never writes files. It answers the
questions you have before pre-processing: what is the crystal's space group,
how many independent elastic constants does it have, which named workflow is
meant for it, and does that workflow work with the crystal in its current
orientation (otherwise: what does).
"""
from dataclasses import dataclass, field

from elastic_vasp import auto, orientation, poscar
from elastic_vasp.systems import SYSTEMS

# (last space group, crystal system, Laue class, named workflow)
_RANGES = [
    (2, "triclinic", "-1", "triclinic"),
    (15, "monoclinic", "2/m", "monoclinic"),
    (74, "orthorhombic", "mmm", "orthogonal"),
    (88, "tetragonal II", "4/m", "tetragonal_2"),
    (142, "tetragonal I", "4/mmm", "tetragonal_1"),
    (148, "trigonal II", "-3", "trigonal_2"),
    (167, "trigonal I", "-3m", "trigonal_1"),
    (194, "hexagonal", "6/mmm", "hexagonal"),
    (230, "cubic", "m-3m", "cubic"),
]
# trigonal I groups whose 2-fold axes are perpendicular to a in the standard
# setting ("31m" type): C15 instead of C14 appears, so trigonal_2 is the named choice
TRIGONAL_31M = frozenset({149, 151, 153, 157, 159, 162, 163})
# named workflows whose orientation check writes POSCAR.reoriented
REORIENTS = frozenset({"hexagonal", "trigonal_1", "trigonal_2", "monoclinic"})


def classify(number):
    """(crystal system, Laue class, named workflow) for a space-group number 1-230."""
    if not 1 <= number <= 230:
        raise ValueError("space group number must be 1-230, got {}".format(number))
    for last, system, laue, named in _RANGES:
        if number <= last:
            if named == "trigonal_1" and number in TRIGONAL_31M:
                named = "trigonal_2"
            return system, laue, named


@dataclass
class SymmetryReport:
    """Result of check_symmetry()."""

    poscar: str
    formula: str
    natoms: int
    space_group: str
    number: int
    crystal_system: str
    laue_class: str
    constants: tuple
    recommended: str
    recommended_works: bool
    working: list = field(default_factory=list)
    symprec: float = auto.SYMPREC

    @property
    def n_constants(self):
        return len(self.constants)

    def summary(self):
        lines = [
            "{}: {} ({} atoms)".format(self.poscar, self.formula, self.natoms),
            "Space group:        {} – {} (Laue class {}), spglib symprec = {} A".format(
                self.space_group, self.crystal_system, self.laue_class, self.symprec),
            "Elastic constants:  {} independent ({} in this orientation)".format(
                self.n_constants, ", ".join(self.constants)),
        ]
        named = "{0}_pre() / {0}_post()".format(self.recommended)
        if self.recommended_works:
            lines.append("Recommended:        {}  – works in this orientation".format(named))
            others = [w for w in self.working if w != self.recommended]
            if others:
                lines.append("Also works:         {} (more calculations)".format(", ".join(others)))
        else:
            lines.append("Recommended:        {}  – but NOT in this orientation: the crystal is "
                         "not aligned the way its deformations expect".format(named))
            if self.working:
                lines.append("Works as it is:     {}".format(", ".join(
                    "{}_pre()".format(w) for w in self.working)))
            if self.recommended in REORIENTS:
                hint = "its warning writes a rotated POSCAR.reoriented"
            else:
                hint = "see the required orientation in the documentation"
            lines.append("To use {}_pre(), rotate the structure first ({}), or use the "
                         "automatic mode.".format(self.recommended, hint))
        lines.append("Any orientation:    auto_pre() / auto_post() – {} deformations".format(
            self.n_constants))
        return "\n".join(lines)


def check_symmetry(poscar_path="POSCAR", symprec=auto.SYMPREC, verbose=True):
    """Report the crystal symmetry of a POSCAR and which workflow to use.

    Prints a short summary (unless ``verbose=False``) and returns a
    :class:`SymmetryReport`. Nothing is written to disk.
    """
    import spglib

    atoms = poscar.read(poscar_path)[0]
    cell = (atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers)
    dataset = spglib.get_symmetry_dataset(cell, symprec=symprec)
    if dataset is None:
        raise ValueError("spglib could not determine the symmetry of {} (symprec={})".format(
            poscar_path, symprec))
    number = dataset.number if hasattr(dataset, "number") else dataset["number"]
    system, laue, named = classify(number)
    basis = auto.invariant_basis(auto.crystal_rotations(atoms, symprec))
    working = sorted((name for name, table in SYSTEMS.items()
                      if orientation.table_covers(table, basis)),
                     key=lambda name: (len(SYSTEMS[name].constants), name))
    report = SymmetryReport(
        poscar=str(poscar_path), formula=atoms.get_chemical_formula(), natoms=len(atoms),
        space_group=orientation.space_group(atoms, symprec), number=int(number),
        crystal_system=system, laue_class=laue, constants=tuple(basis),
        recommended=named, recommended_works=named in working, working=working,
        symprec=symprec)
    if verbose:
        print(report.summary())
    return report
