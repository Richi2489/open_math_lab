"""
metrics.py — métricas de espaciados y correlación de pares
==========================================================

- Distribución de espaciados de vecino más cercano (histograma + CDF).
- Correlación de pares R_2 (la firma Montgomery–Odlyzko) vs 1 − (sin(πx)/(πx))².
- Distancias KS y Wasserstein entre distribuciones de espaciados.
- Momentos de los espaciados.
- Estabilidad por ventana de altura (¿el match GUE es estable o hay efectos de tamaño
  finito a baja altura?).

KS y Wasserstein vienen de scipy.stats (no se reimplementan).
"""

from __future__ import annotations

import numpy as np
from scipy.stats import ks_2samp, wasserstein_distance

from . import gue, spacing


# ---------------------------------------------------------------------------
# Distribución de espaciados
# ---------------------------------------------------------------------------
def histograma(spac, bins=60, rango=(0.0, 4.0)):
    """Histograma normalizado (densidad) de espaciados. Devuelve (centros, densidad)."""
    dens, bordes = np.histogram(spac, bins=bins, range=rango, density=True)
    centros = 0.5 * (bordes[:-1] + bordes[1:])
    return centros, dens


def cdf_empirica(spac):
    """CDF empírica. Devuelve (x ordenado, F)."""
    x = np.sort(np.asarray(spac, dtype=np.float64))
    F = np.arange(1, len(x) + 1) / len(x)
    return x, F


def momentos(spac) -> dict:
    """Media, varianza, asimetría y P(s<0.5) (indicador de repulsión)."""
    s = np.asarray(spac, dtype=np.float64)
    media = float(s.mean())
    var = float(s.var())
    sk = float(((s - media) ** 3).mean() / var**1.5) if var > 0 else 0.0
    return {"media": media, "var": var, "asimetria": sk,
            "frac_menor_0.5": float((s < 0.5).mean())}


# ---------------------------------------------------------------------------
# Distancias entre distribuciones
# ---------------------------------------------------------------------------
def distancias(spac_a, spac_b) -> dict:
    """KS (estadístico y p) y Wasserstein entre dos muestras de espaciados."""
    ks = ks_2samp(spac_a, spac_b)
    return {"ks": float(ks.statistic), "ks_p": float(ks.pvalue),
            "wasserstein": float(wasserstein_distance(spac_a, spac_b))}


# ---------------------------------------------------------------------------
# Correlación de pares R_2 (firma Montgomery–Odlyzko)
# ---------------------------------------------------------------------------
def r2_teorico_gue(r):
    """R_2(r) = 1 − (sin(πr)/(πr))² para GUE."""
    r = np.asarray(r, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        sinc = np.where(r == 0, 1.0, np.sin(np.pi * r) / (np.pi * r))
    return 1.0 - sinc**2


def correlacion_pares(unfolded, r_max=3.0, dr=0.1):
    """Correlación de pares empírica de puntos unfolded (densidad media 1).

    El conteo de pares ordenados (i<j) con separación ≤ r se obtiene VECTORIZADO con
    búsqueda binaria: C(r) = Σ_i (#{j>i : w_j ≤ w_i + r}). El histograma en [r, r+dr] es
    C(r+dr) − C(r), normalizado por (N·dr) = conteo esperado bajo densidad 1 sin
    correlación. Escala a millones de puntos (≈ n_bins búsquedas binarias sobre N).
    Devuelve (centros_r, R_2_empírico).
    """
    w = np.sort(np.asarray(unfolded, dtype=np.float64))
    N = len(w)
    bordes = np.arange(0.0, r_max + dr, dr)
    indices = np.arange(N)

    def C(r):
        hi = np.searchsorted(w, w + r, side="right")
        return float(np.clip(hi - indices - 1, 0, None).sum())

    acum = np.array([C(r) for r in bordes])
    conteo = np.diff(acum)
    centros = 0.5 * (bordes[:-1] + bordes[1:])
    return centros, conteo / (N * dr)


# ---------------------------------------------------------------------------
# Estabilidad por ventana de altura
# ---------------------------------------------------------------------------
def estabilidad_por_altura(gammas, n_ventanas, rng, n_gue=400, recorte=0.1) -> list:
    """Parte los ceros en `n_ventanas` bloques contiguos por altura y, en cada uno,
    compara los espaciados unfolded contra GUE finito simulado (mismo nº de espaciados).

    Devuelve una lista de dicts por ventana: rango de T, n, KS y Wasserstein vs GUE.
    """
    gammas = np.asarray(gammas, dtype=np.float64)
    bloques = np.array_split(np.arange(len(gammas)), n_ventanas)
    res = []
    for idx in bloques:
        if len(idx) < 20:
            continue
        g = gammas[idx]
        s_zeta = spacing.gaps_vecino(g)
        # GUE con ~el mismo número de espaciados que esta ventana.
        n_mat = max(1, int(np.ceil(len(s_zeta) / (n_gue * (1 - 2 * recorte)))))
        s_gue = gue.espaciados_gue(n_gue, rng, n_matrices=n_mat, recorte=recorte)
        d = distancias(s_zeta, s_gue)
        res.append({"T_min": float(g[0]), "T_max": float(g[-1]),
                    "n_spac": int(len(s_zeta)), **d})
    return res
