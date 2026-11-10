"""
Created on Sat Feb 29 16:07:35 2020

@author: pranav
"""
import numpy as np

from elastic_vasp import poscar


def calc_gen(d1, d2, d3, d4, d5, d6, poscar_path="POSCAR", output_path="output.txt"):
    """Write `poscar_path` deformed by the Voigt strain (d1..d6) to `output_path`."""
    lines = poscar.read(poscar_path)

    if poscar.is_plain_scale(lines):
        # a single positive scale factor commutes with the strain: keep the file as written
        lattice = poscar.raw_lattice(lines)
        scale_line = lines[1]
    else:
        # A negative scale (target volume) or per-axis scales would make VASP rescale
        # the strained cell and undo part of the strain, so write the actual
        # lattice with scale 1.0 instead.
        lattice = poscar.lattice(lines)
        scale_line = "1.0\n"

    strain_matrix = np.array(
        [[1 + d1, d6 / 2, d5 / 2],
         [d6 / 2, 1 + d2, d4 / 2],
         [d5 / 2, d4 / 2, 1 + d3]]
    )
    out = np.matmul(lattice, strain_matrix)

    with open(output_path, "w") as my_file:
        my_file.write(lines[0])
        my_file.write(scale_line)
        np.savetxt(my_file, out, fmt="%10.12f")
        my_file.writelines(lines[5:])
