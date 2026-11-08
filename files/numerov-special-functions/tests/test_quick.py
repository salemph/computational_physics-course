"""Fast regression tests: each checks one algorithm against an exact result.
Run with `pytest -q tests` (or `make test`)."""
import importlib.util, pathlib, sys
import numpy as np
from scipy.special import jv, spherical_jn, eval_genlaguerre

ROOT = pathlib.Path(__file__).resolve().parent.parent


def load(rel):
    path = ROOT / rel
    spec = importlib.util.spec_from_file_location(path.stem, path)
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(path.parent))
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- volume 1
def test_free_particle_box_E1():
    p1 = load("vol1_numerov/python/p1_free_particle.py")
    E, Eex, *_ = p1.box_eigen(1, N=200)
    assert abs(E - Eex) / Eex < 1e-8


def test_numerov_fourth_order():
    p1 = load("vol1_numerov/python/p1_free_particle.py")
    e1 = p1.part_A(N=100)[2]; e2 = p1.part_A(N=200)[2]
    assert 14 < e1 / e2 < 18                     # ratio 16 <=> order 4


def test_harmonic_oscillator_levels():
    p3 = load("vol1_numerov/python/p3_harmonic_oscillator.py")
    _, _, st = p3.eigenstates(nmax=2, L=7.0, N=700, Emax=3.0)
    for n, E, _ in st:
        assert abs(E - (n + 0.5)) < 1e-7


def test_square_barrier_transmission():
    p2 = load("vol1_numerov/python/p2_square_barrier.py")
    T, R = p2.scatter(2.0, 0.01)
    assert abs(T - p2.T_exact(2.0)) < 1e-8 and abs(R + T - 1) < 1e-8


# ---------------------------------------------------------------- volume 2
def test_laguerre_recurrence():
    c1 = load("vol2_special_functions/python/ch1_laguerre_hydrogen.py")
    x = np.linspace(0.1, 30, 50)
    L = c1.laguerre_rec(12, 2.0, x)
    assert np.max(np.abs(L[12] - eval_genlaguerre(12, 2.0, x)) / (1 + np.abs(L[12]))) < 1e-12


def test_hydrogen_2p():
    c1 = load("vol2_special_functions/python/ch1_laguerre_hydrogen.py")
    E, u, r = c1.eigen(2, 1, 0.02)
    assert abs(E + 0.125) < 1e-9


def test_legendre_orthonormality():
    c2 = load("vol2_special_functions/python/ch2_legendre_rotor.py")
    x, w = c2.gauss_legendre(30)
    P = c2.plm_bar(15, 2, x)[2:]
    assert np.max(np.abs((P * w) @ P.T - np.eye(len(P)))) < 1e-13


def test_angular_eigenvalue():
    c2 = load("vol2_special_functions/python/ch2_legendre_rotor.py")
    lam = c2.eigen(2, 1, 200)[0]
    assert abs(lam - 6) < 1e-7


def test_bessel_miller():
    c3 = load("vol2_special_functions/python/ch3_bessel_circular.py")
    for x in (0.5, 7.3, 25.0):
        assert np.max(np.abs(c3.J_miller(10, x) - jv(np.arange(11), x))) < 1e-14


def test_spherical_bessel_miller():
    c4 = load("vol2_special_functions/python/ch4_spherical_bessel_well.py")
    for x in (0.5, 3.0, 20.0):
        assert np.max(np.abs(c4.j_miller(15, x) - spherical_jn(np.arange(16), x))) < 1e-13
