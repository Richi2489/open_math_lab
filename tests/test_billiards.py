"""
test_billiards.py — pruebas del par integrable/caótico
======================================================

- el solver FD reproduce los autovalores analíticos del rectángulo;
- el conteo cumple la ley de Weyl;
- el rectángulo da estadística ~Poisson y el estadio ~GOE (dentro de banda).
"""

import numpy as np
import pytest

from billiards import solver, weyl
from riemann import gue, metrics


def _espaciados_norm(E, area, perim, descartar=0):
    s = weyl.espaciados(np.sort(np.asarray(E, float))[descartar:], area, perim)
    return s / s.mean()


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
def test_rectangulo_fd_reproduce_analitico():
    a, b = 1.0, 1.3
    ana = solver.eigs_rectangulo(a, b, 20)
    fd = solver.resolver_dirichlet(
        lambda X, Y: solver.rectangulo_dentro(X, Y, a, b), 0, a, 0, b, 0.012, 20)
    err = np.abs(fd[:15] - ana[:15]) / ana[:15]
    assert err.max() < 0.02  # < 2% en los modos bajos


def test_eigs_rectangulo_ordenados_y_formula():
    a, b = 1.0, np.sqrt(2)
    E = solver.eigs_rectangulo(a, b, 50)
    assert np.all(np.diff(E) >= 0)
    e11 = np.pi ** 2 * (1 / a ** 2 + 1 / b ** 2)  # modo fundamental (1,1)
    assert abs(E[0] - e11) < 1e-9


# ---------------------------------------------------------------------------
# Ley de Weyl
# ---------------------------------------------------------------------------
def test_weyl_conteo():
    a, b = 1.0, np.sqrt(2)
    E = solver.eigs_rectangulo(a, b, 1000)
    area, perim = weyl.geom_rectangulo(a, b)
    # N_weyl(E_1000) debe estar cerca de 1000 (área + perímetro)
    assert abs(weyl.N_weyl(E[-1], area, perim) - 1000) / 1000 < 0.05


# ---------------------------------------------------------------------------
# Estadística: rectángulo ~ Poisson, estadio ~ GOE
# ---------------------------------------------------------------------------
_GOLDEN = (1.0 + np.sqrt(5)) / 2.0


def test_rectangulo_es_poisson():
    # razón áurea: cuadrado irracional, evita las degeneraciones aritméticas de √2.
    a, b = 1.0, _GOLDEN
    E = solver.eigs_rectangulo(a, b, 2000)
    area, perim = weyl.geom_rectangulo(a, b)
    s = _espaciados_norm(E, area, perim, descartar=30)
    assert abs(s.var() - 1.0) < 0.2          # Poisson: var ≈ 1
    assert (s < 0.1).mean() > 0.07           # sin hoyo de repulsión


def test_goe_baseline_var():
    rng = np.random.default_rng(0)
    s = gue.espaciados_goe(200, rng, n_matrices=20)
    assert abs(s.var() - 0.286) < 0.05       # GOE teórico ≈ 0.286


def test_estadio_es_goe():
    # estadio grueso y pequeño (rápido) pero suficiente para ver repulsión.
    a, r = 1.0, 1.0
    E = solver.resolver_dirichlet(
        lambda X, Y: solver.cuarto_estadio_dentro(X, Y, a, r), 0, a + r, 0, r, 0.03, 200)
    area, perim = weyl.geom_cuarto_estadio(a, r)
    s = _espaciados_norm(E, area, perim, descartar=20)
    # repulsión: varianza muy por debajo de Poisson (1.0), del orden GOE.
    assert s.var() < 0.5
    assert (s < 0.1).mean() < 0.06           # hoyo de repulsión en s→0

    # y comparado con el rectángulo: el estadio reprime mucho más.
    Er = solver.eigs_rectangulo(a, _GOLDEN, 2000)
    sr = _espaciados_norm(Er, *weyl.geom_rectangulo(a, _GOLDEN), descartar=20)
    assert s.var() < sr.var()
