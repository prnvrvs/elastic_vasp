"""
Created on Sat Jun  5 19:16:26 2021

@author: pranav
"""
import os
import shutil

from elastic_vasp.target import calc_gen


def _read_strains(path="strain.dat"):
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            return [item.strip() for item in line.split(",") if item.strip()]
    raise ValueError("strain.dat does not contain any strain values")


def _normalize_strain_label(raw_strain):
    strain_value = float(raw_strain)
    if abs(strain_value) < 1e-12:
        return "0.00"
    return raw_strain.strip()


def _reset_folder(folder):
    if os.path.isdir(folder):
        shutil.rmtree(folder)
    elif os.path.exists(folder):
        os.remove(folder)


def _square_ratio(value, denominator):
    if abs(denominator) < 1e-15:
        raise ValueError("Invalid strain value {} produces division by zero".format(value))
    return (value ** 2) / denominator


def _generate_preprocessing(folder_specs):
    strains = _read_strains()
    cwd = os.getcwd()

    for folder, deformation in folder_specs:
        _reset_folder(folder)
        os.mkdir(folder)

        for raw_strain in strains:
            strain_value = float(raw_strain)
            strain_label = _normalize_strain_label(raw_strain)
            strain_dir = os.path.join(folder, "strain_" + strain_label)

            os.mkdir(strain_dir)
            calc_gen(*deformation(strain_value))

            output_path = os.path.join(cwd, "output.txt")
            if not os.path.exists(output_path):
                raise IOError("calc_gen did not create output.txt")
            shutil.copy(output_path, os.path.join(strain_dir, "POSCAR"))
            os.remove(output_path)


def cubic_pre():
    _generate_preprocessing([
        ("C11_C12_I", lambda s: (s, -s, 0.0, 0.0, 0.0, 0.0)),
        ("C11_C12_II", lambda s: (s, s, s, 0.0, 0.0, 0.0)),
        ("C44", lambda s: (0.0, 0.0, 0.0, s, s, s)),
    ])
    return


def orthogonal_pre():
    _generate_preprocessing([
        ("C11_C22_C12_I", lambda s: (s, -s, _square_ratio(s, 1 - s ** 2), 0.0, 0.0, 0.0)),
        ("C11_C22_C12_II", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
        ("C22_C33_C23_I", lambda s: (_square_ratio(s, 1 - s ** 2), s, -s, 0.0, 0.0, 0.0)),
        ("C22_C33_C23_II", lambda s: (0.0, s, s, 0.0, 0.0, 0.0)),
        ("C11_C33_C13_I", lambda s: (s, _square_ratio(s, 1 - s ** 2), -s, 0.0, 0.0, 0.0)),
        ("C11_C33_C13_II", lambda s: (s, 0.0, s, 0.0, 0.0, 0.0)),
        ("C44", lambda s: (_square_ratio(s, 4 - s ** 2), 0.0, 0.0, s, 0.0, 0.0)),
        ("C55", lambda s: (0.0, _square_ratio(s, 4 - s ** 2), 0.0, 0.0, s, 0.0)),
        ("C66", lambda s: (0.0, 0.0, _square_ratio(s, 4 - s ** 2), 0.0, 0.0, s)),
    ])
    return


def hexagonal_pre():
    _generate_preprocessing([
        ("C11_C12_I", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
        ("C11_C12_II", lambda s: (s, -s, _square_ratio(s, 1 - s ** 2), 0.0, 0.0, 0.0)),
        ("C11_C33_C13_I", lambda s: (s, _square_ratio(s, 1 - s ** 2), -s, 0.0, 0.0, 0.0)),
        ("C44", lambda s: (0.0, 0.0, 0.0, s, s, 0.0)),
        ("C11_C12_C13_C33", lambda s: (s, s, s, 0.0, 0.0, 0.0)),
    ])
    return


def trigonal_1_pre():
    _generate_preprocessing([
        ("C11_C12_I", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
        ("C11_C12_II", lambda s: (s, -s, _square_ratio(s, 1 - s ** 2), 0.0, 0.0, 0.0)),
        ("C11_C33_C13_I", lambda s: (s, 0.0, s, 0.0, 0.0, 0.0)),
        ("C11_C33_C13_II", lambda s: (s, _square_ratio(s, 1 - s ** 2), -s, 0.0, 0.0, 0.0)),
        ("C11_C44_C14_I", lambda s: (s, 0.0, 0.0, s, 0.0, 0.0)),
        ("C11_C44_C14_II", lambda s: (s, 0.0, 0.0, -s, 0.0, 0.0)),
    ])
    return


def trigonal_2_pre():
    _generate_preprocessing([
        ("C11_C12_I", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
        ("C11_C12_II", lambda s: (s, -s, _square_ratio(s, 1 - s ** 2), 0.0, 0.0, 0.0)),
        ("C11_C33_C13_I", lambda s: (s, 0.0, s, 0.0, 0.0, 0.0)),
        ("C11_C33_C13_II", lambda s: (s, _square_ratio(s, 1 - s ** 2), -s, 0.0, 0.0, 0.0)),
        ("C11_C44_C14_I", lambda s: (s, 0.0, 0.0, s, 0.0, 0.0)),
        ("C11_C44_C14_II", lambda s: (s, 0.0, 0.0, -s, 0.0, 0.0)),
        ("C11_C44_C15", lambda s: (s, 0.0, 0.0, 0.0, -s, 0.0)),
    ])
    return


def monoclinic_pre():
    _generate_preprocessing([
        ("C11_C22_C12_I", lambda s: (s, -s, _square_ratio(s, 1 - s ** 2), 0.0, 0.0, 0.0)),
        ("C11_C22_C12_II", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
        ("C22_C33_C23_I", lambda s: (_square_ratio(s, 1 - s ** 2), s, -s, 0.0, 0.0, 0.0)),
        ("C22_C33_C23_II", lambda s: (0.0, s, s, 0.0, 0.0, 0.0)),
        ("C11_C33_C13_I", lambda s: (s, _square_ratio(s, 1 - s ** 2), -s, 0.0, 0.0, 0.0)),
        ("C11_C33_C13_II", lambda s: (s, 0.0, s, 0.0, 0.0, 0.0)),
        ("C44", lambda s: (_square_ratio(s, 4 - s ** 2), 0.0, 0.0, s, 0.0, 0.0)),
        ("C11_C55_C15_I", lambda s: (s, 0.0, 0.0, 0.0, s, 0.0)),
        ("C11_C55_C15_II", lambda s: (s, 0.0, 0.0, 0.0, -s, 0.0)),
        ("C22_C55_C25", lambda s: (0.0, s, 0.0, 0.0, s, 0.0)),
        ("C33_C55_C35", lambda s: (0.0, 0.0, s, 0.0, s, 0.0)),
        ("C44_C66_C46_I", lambda s: (0.0, 0.0, 0.0, s, 0.0, s)),
        ("C44_C66_C46_II", lambda s: (0.0, 0.0, 0.0, 0.0, s, -s)),
    ])
    return


def triclinic_pre():
    _generate_preprocessing([
        ("C11", lambda s: (s, 0.0, 0.0, 0.0, 0.0, 0.0)),
        ("C11_C22_C12", lambda s: (s, s, 0.0, 0.0, 0.0, 0.0)),
        ("C22", lambda s: (0.0, s, 0.0, 0.0, 0.0, 0.0)),
        ("C22_C33_C23", lambda s: (0.0, s, s, 0.0, 0.0, 0.0)),
        ("C33", lambda s: (0.0, 0.0, s, 0.0, 0.0, 0.0)),
        ("C11_C33_C13", lambda s: (s, 0.0, s, 0.0, 0.0, 0.0)),
        ("C44", lambda s: (0.0, 0.0, 0.0, s, 0.0, 0.0)),
        ("C55", lambda s: (0.0, 0.0, 0.0, 0.0, s, 0.0)),
        ("C66", lambda s: (0.0, 0.0, 0.0, 0.0, 0.0, s)),
        ("C11_C44_C14", lambda s: (s, 0.0, 0.0, s, 0.0, 0.0)),
        ("C11_C55_C15", lambda s: (s, 0.0, 0.0, 0.0, s, 0.0)),
        ("C11_C66_C16", lambda s: (s, 0.0, 0.0, 0.0, 0.0, s)),
        ("C22_C44_C24", lambda s: (0.0, s, 0.0, s, 0.0, 0.0)),
        ("C22_C55_C25", lambda s: (0.0, s, 0.0, 0.0, s, 0.0)),
        ("C22_C66_C26", lambda s: (0.0, s, 0.0, 0.0, 0.0, s)),
        ("C33_C44_C34", lambda s: (0.0, 0.0, s, s, 0.0, 0.0)),
        ("C33_C55_C35", lambda s: (0.0, 0.0, s, 0.0, s, 0.0)),
        ("C33_C66_C36", lambda s: (0.0, 0.0, s, 0.0, 0.0, s)),
        ("C44_C55_C45", lambda s: (0.0, 0.0, 0.0, s, s, 0.0)),
        ("C44_C66_C46", lambda s: (0.0, 0.0, 0.0, s, 0.0, s)),
        ("C55_C66_C56", lambda s: (0.0, 0.0, 0.0, 0.0, s, s)),
    ])
    return
