"""
sampling.py — muestreo de trayectorias del mapa acelerado
=========================================================

Genera secuencias de v's para los experimentos de la iteración 2.

CRÍTICO: el escalamiento por bits arranca en números de hasta 64 bits, cuyos picos
de trayectoria desbordarían numpy int64 EN SILENCIO. Por eso aquí los arranques se
generan con `random.Random` (enteros nativos de Python, precisión arbitraria) y las
trayectorias se calculan con `engine.secuencia_v`, que también usa enteros de Python.
"""

from __future__ import annotations

import random

import numpy as np

from .engine import secuencia_v

# Buckets de longitud L para la estratificación (etiqueta, lo_inclusive, hi_exclusive).
BUCKETS_L = [
    ("L<25", 0, 25),
    ("25-50", 25, 50),
    ("50-100", 50, 100),
    ("100-200", 100, 200),
    ("L>=200", 200, np.inf),
]

# Cinturón de seguridad: ninguna trayectoria de arranque < 2^64 debería acercarse.
_MAX_ITER = 100_000


def semilla_impar_de_bits(b: int, rng: random.Random) -> int:
    """Entero impar con EXACTAMENTE `b` bits (bit más significativo en 1)."""
    if b < 2:
        raise ValueError("b debe ser >= 2.")
    n = (1 << (b - 1)) | rng.getrandbits(b - 1)
    return n | 1  # forzar impar


def muestrear_banda(
    banda: tuple[int, int],
    n_eff_objetivo: int,
    rng: random.Random,
    min_L: int = 11,
):
    """Muestrea trayectorias con arranque en [lo, hi] bits hasta alcanzar el N_eff
    objetivo en lag-1 (es decir, sum_i (L_i - 1) >= n_eff_objetivo).

    Fijar un N_eff objetivo (en vez de un número de trayectorias) hace que las
    bandas de significancia sean comparables ENTRE bandas de bits.

    Solo se conservan trayectorias con L >= min_L (suficientes para los lags pedidos).
    Devuelve (secuencias, longitudes, n_eff_alcanzado, n_descartadas_por_corto).
    """
    lo, hi = banda
    secuencias: list[np.ndarray] = []
    longitudes: list[int] = []
    n_eff = 0
    descartadas = 0
    while n_eff < n_eff_objetivo:
        b = rng.randint(lo, hi)
        n = semilla_impar_de_bits(b, rng)
        vs = secuencia_v(n, max_iter=_MAX_ITER)
        L = len(vs)
        if L < min_L:
            descartadas += 1
            continue
        secuencias.append(vs)
        longitudes.append(L)
        n_eff += L - 1
    return secuencias, np.array(longitudes, dtype=np.int64), n_eff, descartadas


def bucket_de_longitud(L: int) -> str:
    """Etiqueta del bucket de longitud al que pertenece L."""
    for etiqueta, lo, hi in BUCKETS_L:
        if lo <= L < hi:
            return etiqueta
    return BUCKETS_L[-1][0]


def estratificar_por_longitud(secuencias) -> dict:
    """Agrupa las secuencias por bucket de longitud. Devuelve {etiqueta: [secuencias]}."""
    grupos: dict[str, list] = {etiqueta: [] for etiqueta, _, _ in BUCKETS_L}
    for v in secuencias:
        grupos[bucket_de_longitud(len(v))].append(v)
    return grupos
