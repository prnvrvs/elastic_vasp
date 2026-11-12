"""Command line interface: ``elastic-vasp check | pre | post | systems``.

Examples::

    elastic-vasp check                       # space group, constants, which workflow
    elastic-vasp pre cubic                    # strains from strain.dat
    elastic-vasp pre cubic --strains -0.02 -0.01 0 0.01 0.02
    elastic-vasp pre auto --max-strain 0.03 --points 7
    elastic-vasp post cubic --output C.txt    # also save the 6x6 matrix
    elastic-vasp post auto --frame standard -o C.txt   # C in the standard setting
    elastic-vasp -C path/to/calc post auto

Also available as ``python -m elastic_vasp``.
"""
import argparse
import os
import sys
import warnings

SYSTEM_CHOICES = ["auto", "cubic", "hexagonal", "tetragonal_1", "tetragonal_2", "trigonal_1",
                  "trigonal_2", "orthogonal", "monoclinic", "triclinic"]


def _parser():
    from elastic_vasp import __version__
    from elastic_vasp.auto import SYMPREC

    parser = argparse.ArgumentParser(
        prog="elastic-vasp",
        description="Second-order elastic constants from VASP energy-strain calculations.",
        epilog="Documentation: https://prnvrvs.github.io/elastic_vasp/")
    parser.add_argument("--version", action="version", version="%(prog)s " + __version__)
    parser.add_argument("-C", "--directory", default=".",
                        help="run in this folder (default: current folder)")
    sub = parser.add_subparsers(dest="command", required=True, metavar="COMMAND")

    check = sub.add_parser("check", help="report the crystal symmetry and which workflow to use")
    check.add_argument("poscar", nargs="?", default="POSCAR", help="structure file (default: POSCAR)")
    check.add_argument("--symprec", type=float, default=SYMPREC,
                       help="spglib symmetry tolerance in Angstrom (default: %(default)s)")

    pre = sub.add_parser("pre", help="write the strained POSCARs (pre-processing)")
    pre.add_argument("system", choices=SYSTEM_CHOICES, metavar="SYSTEM",
                     help="crystal system or 'auto': " + ", ".join(SYSTEM_CHOICES))
    strains = pre.add_mutually_exclusive_group()
    strains.add_argument("--strains", nargs="+", type=float, metavar="S",
                         help="strain values, e.g. --strains -0.02 -0.01 0 0.01 0.02 "
                              "(default: read strain.dat)")
    strains.add_argument("--max-strain", type=float, metavar="S",
                         help="symmetric strains from -S to +S (use with --points)")
    pre.add_argument("--points", type=int, default=7,
                     help="number of strains with --max-strain (default: %(default)s)")
    pre.add_argument("--overwrite", action="store_true",
                     help="replace existing deformation folders (deletes their contents)")
    pre.add_argument("--symprec", type=float, default=SYMPREC,
                     help="spglib tolerance for 'auto' (default: %(default)s)")

    post = sub.add_parser("post", help="fit the VASP energies and print the elastic constants")
    post.add_argument("system", choices=SYSTEM_CHOICES, metavar="SYSTEM",
                      help="the same SYSTEM as used for 'pre'")
    post.add_argument("-o", "--output", metavar="FILE",
                      help="also save the 6x6 stiffness matrix (GPa) to FILE")
    post.add_argument("--frame", choices=["input", "standard"], default="input",
                      help="for 'auto': save/return C in the orientation of your POSCAR (input, "
                           "default) or in the crystal's standard setting (standard)")

    sub.add_parser("systems", help="list the crystal systems and their deformations")
    return parser


def _show_warning(message, category, filename, lineno, file=None, line=None):
    print("warning: {}".format(message), file=sys.stderr)


def _check(args):
    from elastic_vasp.symmetry import check_symmetry

    check_symmetry(args.poscar, symprec=args.symprec)


def _pre(args):
    from elastic_vasp import elastic_gen as eg

    strains = args.strains
    if args.max_strain is not None:
        strains = eg.symmetric_strains(args.max_strain, args.points)
    if args.system == "auto":
        eg.auto_pre(overwrite=args.overwrite, symprec=args.symprec, strains=strains)
    else:
        getattr(eg, args.system + "_pre")(overwrite=args.overwrite, strains=strains)
    print("Deformation folders written. Run VASP in every strain_* folder, then: "
          "elastic-vasp post {}".format(args.system))


def _post(args):
    import numpy as np

    from elastic_vasp import constants

    if args.system == "auto":
        c = constants.auto_post(frame=args.frame)
    elif args.frame != "input":
        raise ValueError("--frame standard is only for 'auto' (the named workflows already "
                         "work in the standard orientation)")
    else:
        c = constants.analyse(args.system)
    if args.output:
        # plain 6x6 text (GPa, Voigt) without a header line, so tools such as
        # elastic-vis (https://github.com/prnvrvs/elastic-vis) can read it directly
        np.savetxt(args.output, c, fmt="%12.4f")
        print("Stiffness matrix saved to {}".format(args.output))


def _systems(args):
    from elastic_vasp.systems import SYSTEMS

    for name, spec in SYSTEMS.items():
        print("{:<13} {:>2} constants: {}".format(name, len(spec.constants), ", ".join(spec.constants)))
        for k, d in enumerate(spec.deformations):
            print("    {:<16} {:<54} a = {}".format(d.folder, d.describe(), spec.energy_expression(k)))
    print("auto          constants and deformations from the crystal's symmetry (any orientation)")


def main(argv=None):
    args = _parser().parse_args(argv)
    from elastic_vasp.coefficient import MissingEnergyError

    old_cwd = os.getcwd()
    previous = warnings.showwarning
    warnings.showwarning = _show_warning
    try:
        os.chdir(args.directory)
        {"check": _check, "pre": _pre, "post": _post, "systems": _systems}[args.command](args)
    except (FileNotFoundError, FileExistsError, ValueError, MissingEnergyError, ImportError) as exc:
        print("error: {}".format(exc), file=sys.stderr)
        return 1
    finally:
        warnings.showwarning = previous
        os.chdir(old_cwd)
    return 0


if __name__ == "__main__":
    sys.exit(main())
