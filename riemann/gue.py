"""
gue.py — GUE finito simulado y baseline Poisson
===============================================

CLAVE metodológica: la repulsión de los ceros de ζ se compara contra **GUE de tamaño
finito simulado con el mismo N**, NO contra la fórmula asintótica ideal. Toda desviación
de ζ respecto del GUE asintótico que también aparezca en el GUE finito es artefacto de
tamaño finito, no señal.

- GUE: matriz Hermitiana con entradas gaussianas complejas; sus autovalores, tras
  normalizar por √n, siguen la ley del semicírculo en [−2, 2]. Se hace unfold de los
  autovalores con la CDF del semicírculo para obtener espaciados de media 1.
- Poisson: gaps exponenciales i.i.d. (sin repulsión) — el contraste "nulo sin estructura".
"""

from __future__ import annotations

import numpy as np


def matriz_gue(n: int, rng) -> np.ndarray:
    """Matriz GUE n×n: A complejo gaussiano, H = (A + A†)/2 (Hermitiana).

    Con esta normalización, E|H_ij|²=1 (i≠j) y E H_ii²=1, de modo que los autovalores
    de H/√n tienden a la ley del semicírculo en [−2, 2].
    """
    A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    return (A + A.conj().T) / 2.0


def cdf_semicirculo(x):
    """CDF de la ley del semicírculo en [−2, 2]."""
    x = np.clip(np.asarray(x, dtype=np.float64), -2.0, 2.0)
    return 0.5 + (x * np.sqrt(4.0 - x * x)) / (4.0 * np.pi) + np.arcsin(x / 2.0) / np.pi


def espaciados_gue(n: int, rng, n_matrices: int = 1, recorte: float = 0.1) -> np.ndarray:
    """Espaciados de vecino más cercano de GUE finito (n×n), unfolded por semicírculo.

    Se descartan los bordes (`recorte` por lado) para quedarnos con el bulk, donde la
    ley del semicírculo describe bien la densidad. Devuelve todos los espaciados
    concatenados sobre `n_matrices` realizaciones.
    """
    todos = []
    lo, hi = int(recorte * n), int((1.0 - recorte) * n)
    for _ in range(n_matrices):
        H = matriz_gue(n, rng)
        lam = np.linalg.eigvalsh(H) / np.sqrt(n)
        u = n * cdf_semicirculo(lam)
        u.sort()
        todos.append(np.diff(u[lo:hi]))
    return np.concatenate(todos)


def espaciados_poisson(n_spac: int, rng) -> np.ndarray:
    """Baseline Poisson: gaps exponenciales i.i.d. de media 1 (sin repulsión)."""
    return rng.exponential(1.0, size=int(n_spac))


def wigner_gue(s):
    """Sorpresa de Wigner para GUE: P(s) = (32/π²) s² exp(−4 s²/π) (referencia teórica)."""
    s = np.asarray(s, dtype=np.float64)
    return (32.0 / np.pi**2) * s**2 * np.exp(-4.0 * s**2 / np.pi)
