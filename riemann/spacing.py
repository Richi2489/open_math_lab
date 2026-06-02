"""
spacing.py — unfolding por densidad local y gaps de vecino más cercano
======================================================================

EL PASO CRÍTICO del lab, equivalente al confounder de magnitud en Collatz.

La densidad de ceros de ζ crece con la altura T. El número suave de ceros hasta T es
(Riemann–von Mangoldt, parte suave):

    N(T) ≈ (T/2π)·log(T/2π) − T/2π + 7/8

El *unfolding* mapea cada cero γ_n a w_n = N(γ_n). Así los gaps w_{n+1} − w_n tienen
media 1 LOCALMENTE, sin importar la altura. Los espaciados de vecino más cercano son
s_n = w_{n+1} − w_n.

⚠ ADVERTENCIA (el confounder): normalizar dividiendo por la media GLOBAL de los gaps es
INCORRECTO, porque la densidad cambia con la altura — es el mismo tipo de artefacto de
escala que matamos en Collatz. La función `gaps_normalizacion_global_INCORRECTA` existe
solo para EXHIBIR ese error en el reporte, nunca como método válido.
"""

from __future__ import annotations

import numpy as np

DOS_PI = 2.0 * np.pi


def N_suave(T):
    """Parte suave del contador de ceros de Riemann–von Mangoldt, N(T)."""
    T = np.asarray(T, dtype=np.float64)
    t = T / DOS_PI
    return t * np.log(t) - t + 7.0 / 8.0


def N_suave_inversa(w, t_max=1e7):
    """Inverso numérico de N_suave: dado w, devuelve T con N_suave(T)=w.

    Útil para tests (round-trip) y para generar mallas. Vectorizado por brentq escalar.
    """
    from scipy.optimize import brentq

    def _inv(wi):
        return brentq(lambda T: N_suave(T) - wi, DOS_PI + 1e-6, t_max)

    w = np.atleast_1d(np.asarray(w, dtype=np.float64))
    return np.array([_inv(float(wi)) for wi in w])


def unfold(gammas: np.ndarray) -> np.ndarray:
    """Unfold por densidad local (asintótico): w_n = N_suave(γ_n)."""
    return N_suave(np.asarray(gammas, dtype=np.float64))


def theta_unfold(gammas: np.ndarray) -> np.ndarray:
    """Unfold EXACTO: w_n = θ(γ_n)/π + 1 (Riemann–Siegel, vía mpmath).

    Versión directa (sin caché), pensada para pocos valores / tests. Para muestras
    grandes usar `riemann.zeros.cargar_unfold_exacto`, que cachea.
    """
    import mpmath

    return np.array(
        [float(mpmath.siegeltheta(mpmath.mpf(float(g))) / mpmath.pi + 1)
         for g in np.atleast_1d(np.asarray(gammas, dtype=np.float64))]
    )


def gaps_vecino(gammas: np.ndarray) -> np.ndarray:
    """Espaciados de vecino más cercano de los ceros unfolded (media ≈ 1)."""
    return np.diff(unfold(gammas))


def unfold_densidad_local(gammas: np.ndarray) -> np.ndarray:
    """Unfolding por densidad LOCAL, gap a gap (forma estándar para Odlyzko):

        s_n = (γ_{n+1} − γ_n) · log(γ_n / 2π) / (2π)

    El factor log(γ/2π)/(2π) es la densidad media local de ceros a la altura γ, así que
    s_n tiene media ≈ 1 localmente, sin importar la altura. Equivale (a primer orden) a
    diferenciar N_suave, pero se expresa por gap. Es el método por defecto a alta altura.
    """
    g = np.asarray(gammas, dtype=np.float64)
    gaps = np.diff(g)
    dens = np.log(g[:-1] / DOS_PI) / DOS_PI
    return gaps * dens


def gaps_normalizacion_global_INCORRECTA(gammas: np.ndarray) -> np.ndarray:
    """MÉTODO INCORRECTO (solo para exhibir el confounder en el reporte).

    Divide los gaps crudos por su media global. Como la densidad crece con la altura,
    esto NO produce espaciados estacionarios: mezcla regímenes de densidad distintos y
    distorsiona la distribución. No usar como unfolding válido.
    """
    gaps = np.diff(np.asarray(gammas, dtype=np.float64))
    return gaps / gaps.mean()
