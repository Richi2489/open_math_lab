"""
zeros.py — ceros no triviales de ζ con caché incremental
========================================================

Fuente reproducible: `mpmath.zetazero(n)` da el n-ésimo cero no trivial sobre la línea
crítica, 1/2 + i·γ_n. Aquí solo nos interesa la parte imaginaria γ_n (γ_n > 0).

`mpmath` es LENTO (~0.2 s/cero), así que los γ_n se cachean en `data/` y se calculan de
forma INCREMENTAL: si el caché ya tiene M ceros y se piden N>M, solo se computan los N−M
que faltan y se anexan. Reabrir una corrida no recomputa nada.

Hook extensible: `cargar_tabla_odlyzko` queda como placeholder para ingerir, en una
iteración futura, las tablas de Odlyzko de alta altura T (donde el match GUE es más limpio).
NO se descarga nada todavía.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

_DIR_DATOS = Path(__file__).resolve().parents[1] / "data"
_CACHE = _DIR_DATOS / "zeta_zeros_im.npy"
_META = _DIR_DATOS / "zeta_zeros_meta.json"
_DPS = 25  # dígitos de precisión de mpmath; sobra para γ hasta ~10^4


def _cargar_cache() -> np.ndarray:
    if _CACHE.exists():
        return np.load(_CACHE)
    return np.empty(0, dtype=np.float64)


def cargar_ceros(n: int, dps: int = _DPS, verbose: bool = False) -> np.ndarray:
    """Devuelve las partes imaginarias γ_1..γ_n de los ceros de ζ (np.float64).

    Computa de forma incremental solo los que falten respecto del caché y los persiste.
    """
    if n < 1:
        raise ValueError("n debe ser >= 1.")
    cache = _cargar_cache()
    if len(cache) >= n:
        return cache[:n].copy()

    import mpmath

    mpmath.mp.dps = dps
    nuevos = []
    inicio = len(cache) + 1
    for k in range(inicio, n + 1):
        nuevos.append(float(mpmath.zetazero(k).imag))
        if verbose and (k % 100 == 0 or k == n):
            print(f"    ζ-ceros: {k}/{n}")
    todos = np.concatenate([cache, np.array(nuevos, dtype=np.float64)])

    _DIR_DATOS.mkdir(exist_ok=True)
    np.save(_CACHE, todos)
    _META.write_text(json.dumps({"n": int(len(todos)), "dps": dps}), encoding="utf-8")
    return todos[:n].copy()


def validar_ceros(gammas: np.ndarray) -> None:
    """Sanity de los ceros: positivos y estrictamente crecientes."""
    g = np.asarray(gammas, dtype=np.float64)
    if g.size == 0:
        raise ValueError("lista de ceros vacía.")
    if not np.all(g > 0):
        raise ValueError("hay partes imaginarias no positivas.")
    if not np.all(np.diff(g) > 0):
        raise ValueError("los ceros no son estrictamente crecientes.")


def cargar_tabla_odlyzko(path):  # pragma: no cover - placeholder de iteración futura
    """Hook para ingerir tablas de Odlyzko de alta altura T. No implementado aún:
    en la iteración 1 los ceros vienen solo de mpmath (sin descargas externas)."""
    raise NotImplementedError(
        "Ingesta de tablas de Odlyzko: pendiente para una iteración futura "
        "(alta altura T, donde el match GUE es más limpio)."
    )
