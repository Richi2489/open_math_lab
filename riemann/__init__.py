"""
open_math_lab — riemann
=======================

Laboratorio estadístico reproducible de la conexión Montgomery–Odlyzko: ¿los gaps
normalizados (unfolded) de los ceros no triviales de ζ exhiben repulsión tipo GUE?

NO es proof-hunting de la Hipótesis de Riemann. Misma filosofía anti-autoengaño que el
lab de Collatz (ver docs/riemann_lab_plan.md).

Submódulos:
    zeros   : ceros de ζ vía mpmath, con caché incremental en data/.
    spacing : unfolding por densidad local N(T) y gaps de vecino más cercano.
    gue     : GUE finito simulado y baseline Poisson.
    metrics : espaciados, correlación de pares R_2, KS, Wasserstein, estabilidad por altura.
    viz     : figuras (espaciados, R_2, estabilidad).
"""

from . import gue, metrics, odlyzko, spacing, zeros

__all__ = ["zeros", "spacing", "gue", "metrics", "odlyzko"]
