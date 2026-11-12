"""Consistency checks for the crystal-system tables in systems.py."""
import numpy as np
import pytest

from elastic_vasp.systems import SYSTEMS, CrystalSystem, Deformation

EXPECTED_COUNTS = {"cubic": 3, "hexagonal": 5, "tetragonal_1": 6, "tetragonal_2": 7,
                   "trigonal_1": 6, "trigonal_2": 7, "orthogonal": 9, "monoclinic": 13,
                   "triclinic": 21}


@pytest.mark.parametrize("name", sorted(SYSTEMS))
def test_table_is_consistent(name):
    spec = SYSTEMS[name]
    assert len(spec.constants) == EXPECTED_COUNTS[name]
    assert len({d.folder for d in spec.deformations}) == len(spec.deformations)
    design = spec.design_matrix()
    assert np.linalg.matrix_rank(design) == len(spec.constants)

    rng = np.random.default_rng(1)
    values = rng.normal(size=len(spec.constants))
    c = spec.matrix(values)
    assert c.shape == (6, 6)
    np.testing.assert_allclose(c, c.T)  # stiffness matrices are symmetric
    np.testing.assert_allclose(spec.matrix(np.zeros(len(values))), 0)  # linear, no offset
    # every independent constant appears in the matrix
    for j in range(len(values)):
        assert np.any(spec.matrix(np.eye(len(values))[j]) != 0)


@pytest.mark.parametrize("name", sorted(SYSTEMS))
def test_solve_inverts_the_energy_model(name):
    spec = SYSTEMS[name]
    truth = np.random.default_rng(2).uniform(10, 100, size=len(spec.constants))
    c_true = spec.matrix(truth)
    coefficients = {}
    for d in spec.deformations:
        e = np.array(d.strain, float)
        coefficients[d.folder] = 0.5 * e @ c_true @ e
    values, c = spec.solve(coefficients)
    np.testing.assert_allclose([values[n] for n in spec.constants], truth)
    np.testing.assert_allclose(c, c_true)


def test_underdetermined_system_is_rejected():
    spec = CrystalSystem(
        "broken", SYSTEMS["cubic"].constants, SYSTEMS["cubic"].layout,
        SYSTEMS["cubic"].deformations[:2],
    )
    with pytest.raises(ValueError, match="determine only 2 of 3"):
        spec.solve({d.folder: 1.0 for d in spec.deformations})


def test_deformation_matches_original_strain_formulas():
    d = Deformation("x", (1, -1, 0, 0, 0, 0), (2, 1))
    s = -0.02
    assert d.at(s) == (s, -s, s ** 2 / (1 - s ** 2), 0.0, 0.0, 0.0)
    # zero components are +0.0, never -0.0 (keeps POSCAR text identical)
    assert all(str(v) == "0.0" for v in Deformation("y", (0, 0, 0, 1, 0, 0)).at(-0.01)[:3])


@pytest.mark.parametrize("name", sorted(SYSTEMS))
def test_written_lattice_is_L_times_F(name, tmp_path):
    """Every deformation writes exactly L' = L F with F = I + eps, eps the symmetric strain
    tensor of the Voigt strain (off-diagonal entries are half the engineering shears)."""
    from ase.build import bulk
    from elastic_vasp import poscar
    from elastic_vasp.target import calc_gen

    atoms = bulk("Al", "fcc", a=4.04)
    atoms.set_cell(atoms.cell.array + np.array([[0, 0.1, 0.05], [0.07, 0, 0.02], [0.03, 0.04, 0]]),
                   scale_atoms=True)  # a generic, non-symmetric cell
    src, out = tmp_path / "POSCAR", tmp_path / "out"
    poscar.write(src, atoms, "generic")
    lattice = poscar.lattice(src)
    for d in SYSTEMS[name].deformations:
        for delta in (-0.03, 0.02):
            e1, e2, e3, e4, e5, e6 = d.at(delta)
            f = np.eye(3) + np.array([[e1, e6 / 2, e5 / 2], [e6 / 2, e2, e4 / 2], [e5 / 2, e4 / 2, e3]])
            calc_gen(e1, e2, e3, e4, e5, e6, poscar_path=str(src), output_path=str(out))
            np.testing.assert_allclose(poscar.lattice(out), lattice @ f, atol=1e-12, err_msg=d.folder)
