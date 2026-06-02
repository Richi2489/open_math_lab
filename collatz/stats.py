"""
stats.py — medidas empíricas de la heurística de Collatz
========================================================

Cuatro medidas, alineadas con la pregunta de la iteración 1:

1. E[v] empírico            -> debe ≈ 2.0
   La heurística modela v como una geométrica: P(v=k) = 2^-k, k>=1.
   Entonces E[v] = sum_{k>=1} k·2^-k = 2.

2. Drift geométrico por paso -> debe ≈ 0.75 = 3/4
   Para el mapa acelerado T(n) = (3n+1)/2^v, el factor por paso es T(n)/n ≈ 3/2^v.
   La media GEOMÉTRICA es exp(E[log(T/n)]) = exp(log 3 - E[v]·log 2) = 3/4 si E[v]=2.

3. Autocorrelación de la secuencia de v's de una trayectoria.
   La heurística asume que los v's son i.i.d. (independientes). Aquí lo medimos
   directamente: si hubiera independencia, la autocorrelación sería ~0 para lag>=1.

4. Distribución de tiempos totales de parada.

Todas las funciones aceptan un `rng` (np.random.Generator) para reproducibilidad.
"""

from __future__ import annotations

import numpy as np

from .engine import paso_acelerado_vector, secuencia_v, tiempo_total_parada

LOG3_SOBRE_4 = np.log(0.75)  # valor teórico del drift logarítmico


# ---------------------------------------------------------------------------
# Muestreo
# ---------------------------------------------------------------------------
def muestra_impares(n_muestras: int, bits: int = 40, rng=None) -> np.ndarray:
    """Genera `n_muestras` enteros impares aleatorios de ~`bits` bits (int64).

    Con bits<=40, el cálculo 3n+1 cabe en int64 sin desbordar.
    """
    if bits > 61:
        raise ValueError("bits demasiado grande: 3n+1 desbordaría int64. Usa bits<=61.")
    rng = rng if rng is not None else np.random.default_rng()
    n = rng.integers(1, 1 << bits, size=n_muestras, dtype=np.int64)
    n |= 1  # forzar impar
    return n


# ---------------------------------------------------------------------------
# 1) E[v]
# ---------------------------------------------------------------------------
def esperanza_v(n_muestras: int = 1_000_000, bits: int = 40, rng=None):
    """E[v] empírico sobre un paso acelerado de impares aleatorios.

    Devuelve (media, array_de_v).  Esperado: media ≈ 2.0.
    """
    n = muestra_impares(n_muestras, bits, rng)
    _, v = paso_acelerado_vector(n)
    return float(v.mean()), v


# ---------------------------------------------------------------------------
# 2) Drift geométrico
# ---------------------------------------------------------------------------
def drift_geometrico(n_muestras: int = 1_000_000, bits: int = 40, rng=None):
    """Media geométrica del factor T(n)/n del mapa acelerado.

    Devuelve (media_geometrica, log_ratios).  Esperado: media_geometrica ≈ 0.75.
    """
    n = muestra_impares(n_muestras, bits, rng)
    siguiente, _ = paso_acelerado_vector(n)
    log_ratio = np.log(siguiente.astype(np.float64)) - np.log(n.astype(np.float64))
    media_geom = float(np.exp(log_ratio.mean()))
    return media_geom, log_ratio


# ---------------------------------------------------------------------------
# 3) Autocorrelación de los v's
# ---------------------------------------------------------------------------
def autocorrelacion(serie: np.ndarray, max_lag: int = 40) -> np.ndarray:
    """Función de autocorrelación (ACF) de una serie 1D, lags 0..max_lag.

    Normalizada de modo que acf[0] = 1.
    """
    x = np.asarray(serie, dtype=np.float64)
    x = x - x.mean()
    var = np.dot(x, x)
    if var == 0:
        # serie constante: no hay variabilidad que correlacionar
        acf = np.zeros(max_lag + 1)
        acf[0] = 1.0
        return acf
    max_lag = min(max_lag, len(x) - 1)
    acf = np.empty(max_lag + 1)
    acf[0] = 1.0
    for k in range(1, max_lag + 1):
        acf[k] = np.dot(x[:-k], x[k:]) / var
    return acf


def autocorrelacion_v(n_semilla: int, max_lag: int = 40) -> tuple[np.ndarray, np.ndarray]:
    """ACF de la secuencia de v's de UNA trayectoria acelerada que arranca en `n_semilla`.

    Devuelve (lags, acf).
    """
    vs = secuencia_v(n_semilla)
    acf = autocorrelacion(vs, max_lag=max_lag)
    return np.arange(len(acf)), acf


def autocorrelacion_v_promedio(
    n_trayectorias: int = 2000,
    bits: int = 40,
    max_lag: int = 30,
    min_longitud: int = 40,
    rng=None,
):
    """ACF promediada sobre muchas trayectorias (estimación estable).

    Para cada semilla impar aleatoria se calcula la ACF de sus v's y se promedian
    solo las trayectorias con al menos `min_longitud` pasos (para que el lag tenga
    soporte). Devuelve (lags, acf_promedio, n_usadas).
    """
    rng = rng if rng is not None else np.random.default_rng()
    semillas = muestra_impares(n_trayectorias, bits, rng)
    acumulado = np.zeros(max_lag + 1)
    usadas = 0
    for s in semillas:
        vs = secuencia_v(int(s))
        if len(vs) <= min_longitud:
            continue
        acumulado += autocorrelacion(vs, max_lag=max_lag)
        usadas += 1
    if usadas == 0:
        raise RuntimeError("Ninguna trayectoria alcanzó la longitud mínima requerida.")
    return np.arange(max_lag + 1), acumulado / usadas, usadas


# ---------------------------------------------------------------------------
# 4) Distribución de tiempos de parada
# ---------------------------------------------------------------------------
def distribucion_tiempos_parada(
    n_muestras: int = 20_000, n_max: int = 10**6, rng=None
) -> np.ndarray:
    """Tiempos totales de parada para `n_muestras` semillas aleatorias en [2, n_max)."""
    rng = rng if rng is not None else np.random.default_rng()
    semillas = rng.integers(2, n_max, size=n_muestras, dtype=np.int64)
    return np.array([tiempo_total_parada(int(s)) for s in semillas], dtype=np.int64)
