"""
Created on Sat Feb 29 16:07:35 2020

@author: pranav
"""
import numpy as np

from elastic_vasp import poscar


def calc_gen(d1, d2, d3, d4, d5, d6, poscar_path="POSCAR", output_path="output.txt"):
    """Write `poscar_path` deformed by the Voigt strain (d1..d6) to `output_path`.

    The lattice rows are multiplied by (I + strain) and the atoms keep their
    fractional coordinates, so they move with the cell.
    """
    atoms, title = poscar.read(poscar_path)
    strain_matrix = np.array(
        [[1 + d1, d6 / 2, d5 / 2],
         [d6 / 2, 1 + d2, d4 / 2],
         [d5 / 2, d4 / 2, 1 + d3]]
    )
    atoms.set_cell(atoms.cell.array @ strain_matrix, scale_atoms=True)
    poscar.write(output_path, atoms, title)
