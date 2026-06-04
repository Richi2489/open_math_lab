"""
weyl.py — unfolding por la ley de Weyl (área + perímetro)
=========================================================

La densidad media de modos de un billar 2D Dirichlet sigue la ley de Weyl:

    N(E) ≈ (A/4π)·E − (P/4π)·√E        (área A, perímetro P; término de perímetro negativo
                                        para Dirichlet)

Unfolding correcto = usar AMBOS términos (no solo el de área): w_i = N(E_i). Así los
espaciados w_{i+1}−w_i tienen media ≈ 1 y la estadística es comparable entre dominios.
"""

from __future__ import annotations

import numpy as np

CUATRO_PI = 4.0 * np.pi


def N_weyl(E, area, perim):
    """Conteo suave de Weyl con término de área y de perímetro."""
    E = np.asarray(E, dtype=np.float64)
    return area / CUATRO_PI * E - perim / CUATRO_PI * np.sqrt(E)


def unfold(E, area, perim):
    """Coordenadas unfolded w_i = N_weyl(E_i)."""
    return N_weyl(np.asarray(E, dtype=np.float64), area, perim)


def espaciados(E, area, perim):
    """Espaciados de vecino más cercano tras el unfolding de Weyl (media ≈ 1)."""
    return np.diff(unfold(np.sort(np.asarray(E, dtype=np.float64)), area, perim))


def geom_rectangulo(a, b):
    """(área, perímetro) del rectángulo [0,a]×[0,b]."""
    return a * b, 2.0 * (a + b)


def geom_cuarto_estadio(a, r):
    """(área, perímetro) del cuarto de estadio desimetrizado (Dirichlet en todo el borde,
    ejes de simetría incluidos)."""
    area = a * r + np.pi * r ** 2 / 4.0
    perim = 2.0 * a + 2.0 * r + np.pi * r / 2.0
    return area, perim
