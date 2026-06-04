"""
open_math_lab — billiards
=========================

Billares cuánticos: generar con código la conexión caos → repulsión de niveles que el
programa de Riemann tomó prestada de la teoría de matrices aleatorias
(Bohigas–Giannoni–Schmit, 1984). NO es teoría nueva: es una demostración controlada.

Experimento central (iteración 1): un PAR controlado con el MISMO solver y la MISMA
maquinaria estadística (reusada del lab de Riemann), donde sólo cambia la geometría de la
mesa:

  - rectángulo de lados inconmensurables (integrable)  → estadística de Poisson;
  - cuarto de estadio de Bunimovich (caótico, desimetrizado) → estadística GOE.

Submódulos:
    solver : autovalores del Laplaciano Dirichlet (diferencias finitas) y rectángulo analítico.
    weyl   : unfolding por la ley de Weyl (área + perímetro).
    viz    : figuras (el par de distribuciones de espaciado, CDF, conteo de Weyl).
"""

from . import solver, viz, weyl

__all__ = ["solver", "weyl", "viz"]
