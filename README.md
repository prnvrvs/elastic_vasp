# elastic_vasp

Second-order **elastic constants from VASP** with the energy–strain method. elastic_vasp
writes the strained structures, checks that your crystal and its orientation fit the
deformations, reads the VASP energies back, fits the energy–strain curves and returns the
full 6×6 stiffness matrix with stability checks and Voigt–Reuss–Hill averages.

**Documentation: https://prnvrvs.github.io/elastic_vasp/**

**Video tutorial: [Elastic constants with VASP and elastic_vasp](https://youtu.be/kXaeIoNdHzw)** (8 min, fcc Al)
(legacy video for older versions: [YouTube](https://youtu.be/9z455a8NIME))

## Installation

```bash
pip install "git+https://github.com/prnvrvs/elastic_vasp.git"   # current version (documented)
pip install elastic-vasp                                          # last PyPI release (older)
```

Requires Python 3.10 or newer; NumPy, matplotlib, ASE and spglib are installed
automatically.

## Usage

```python
from elastic_vasp import elastic_gen as eg, constants

eg.cubic_pre()             # 1. write the strained POSCARs
# 2. run VASP (ISIF = 2) in every strain folder
C = constants.cubic_post() # 3. 6x6 stiffness matrix in GPa
```

Or from the shell:

```bash
elastic-vasp check                                 # space group, constants, which workflow
elastic-vasp pre cubic --max-strain 0.03 --points 7
elastic-vasp post cubic --output C.txt
```

Named workflows exist for every crystal system (`cubic`, `hexagonal`, `tetragonal_1`,
`tetragonal_2`, `trigonal_1`, `trigonal_2`, `orthogonal`, `monoclinic`, `triclinic`).
`auto_pre()` / `auto_post()` work out the constants and deformations from the crystal's
symmetry in any orientation. See the
[crystal systems](https://prnvrvs.github.io/elastic_vasp/guide/systems.html) and
[automatic mode](https://prnvrvs.github.io/elastic_vasp/guide/auto.html) pages.

## Visualisation

Plot directional Young's modulus, shear modulus, Poisson's ratio and linear compressibility
with the companion package [elastic-vis](https://github.com/prabhat7iitm/elastic-vis):

```bash
elastic-vasp post cubic --output C.txt
elastic-vis C.txt --plot --Young --Shear --2D --3D
```

## Legacy version

The original code (0.0.4) remains available on the
[`legacy`](https://github.com/prnvrvs/elastic_vasp/tree/legacy) branch and on PyPI
(`pip install "elastic-vasp==0.0.4"`) for reproducing earlier results. It has known errors
for trigonal II, monoclinic and triclinic crystals that are fixed in the current version;
see [Legacy version](https://prnvrvs.github.io/elastic_vasp/legacy.html).

## Citation

If you find this package useful, please cite:

1. *Kumar, P., and Adlakha, I. (July 26, 2022). "Effect of Interstitial Hydrogen on Elastic
   Behavior of Metals: an Ab-Initio Study." ASME. J. Eng. Mater. Technol.
   doi: https://doi.org/10.1115/1.4055097*
2. *Mishra, P., Kumar, P., Neelakantan, L., and Adlakha, I. (2022). "First-principles
   prediction of electrochemical polarization and mechanical behavior in Mg based
   intermetallics." Computational Materials Science 214, 111667.
   doi: https://doi.org/10.1016/j.commatsci.2022.111667*

```
@article{Elastic_vasp,
    author = {Kumar, P. and Adlakha, I.},
    title = "{Effect of Interstitial Hydrogen on Elastic Behavior of Metals: An Ab-Initio Study}",
    journal = {Journal of Engineering Materials and Technology},
    volume = {145},
    number = {1},
    year = {2022},
    month = {08},
    issn = {0094-4289},
    doi = {10.1115/1.4055097},
}

@article{Mishra_2022,
    author = {Mishra, Pragyandipta and Kumar, Pranav and Neelakantan, Lakshman and Adlakha, Ilaksh},
    title = "{First-principles prediction of electrochemical polarization and mechanical behavior in Mg based intermetallics}",
    journal = {Computational Materials Science},
    volume = {214},
    pages = {111667},
    year = {2022},
    month = {11},
    issn = {0927-0256},
    doi = {10.1016/j.commatsci.2022.111667},
}
```

## Contact

https://sites.google.com/view/prnvkmr4
