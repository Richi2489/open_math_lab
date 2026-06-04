"""
explicit_formula.py — reconstruir el conteo de primos desde los ceros de ζ
==========================================================================

Visualización pedagógica clásica (Riemann, 1859): los ceros no triviales de ζ codifican
los primos. NO es un intento de RH.

- ψ(x) (función de Chebyshev / von Mangoldt): ψ(x) = Σ_{p^k ≤ x} ln(p) = Σ_{n≤x} Λ(n).
  Es una escalera que salta ln(p) en cada potencia de primo p^k.

- Fórmula explícita (truncada a N ceros):

    ψ_N(x) = x − 2·√x · Σ_{n=1..N} Re( exp(i·γ_n·ln x) / (½ + i·γ_n) )
                 − ln(2π) − ½·ln(1 − x⁻²)

  El término dominante es `x` (densidad suave de primos, ~PNT); la suma sobre los ceros
  γ_n son las correcciones OSCILANTES que "enganchan" la curva suave a los saltos en los
  primos. Cada cero aporta una frecuencia.

Con N finito aparece el fenómeno de Gibbs (overshoot cerca de cada salto): es el
comportamiento esperado de una suma truncada, NO un error.
"""

from __future__ import annotations

import math

import numpy as np

LOG_2PI = math.log(2.0 * math.pi)


# ---------------------------------------------------------------------------
# Función de von Mangoldt y ψ exacta
# ---------------------------------------------------------------------------
def von_mangoldt(n: int) -> float:
    """Λ(n) = ln(p) si n = p^k (potencia de un primo), 0 en otro caso."""
    n = int(n)
    if n < 2:
        return 0.0
    p = 2
    while p * p <= n:
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return math.log(p) if m == 1 else 0.0
        p += 1
    return math.log(n)  # n es primo


def _criba(n: int) -> np.ndarray:
    """Primos ≤ n (criba de Eratóstenes)."""
    if n < 2:
        return np.empty(0, dtype=np.int64)
    es = np.ones(n + 1, dtype=bool)
    es[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if es[i]:
            es[i * i::i] = False
    return np.nonzero(es)[0].astype(np.int64)


def _saltos_psi(x_max: float):
    """Posiciones (potencias de primo ≤ x_max) ordenadas y su salto ln(p)."""
    pos, salto = [], []
    for p in _criba(int(x_max)):
        lp = math.log(p)
        pk = int(p)
        while pk <= x_max:
            pos.append(pk)
            salto.append(lp)
            pk *= p
    if not pos:
        return np.empty(0), np.empty(0)
    pos = np.array(pos, dtype=np.float64)
    salto = np.array(salto, dtype=np.float64)
    orden = np.argsort(pos)
    return pos[orden], salto[orden]


def psi_exacta(x, x_max=None):
    """ψ(x) exacta (escalera de Chebyshev). Acepta escalar o array."""
    x_arr = np.atleast_1d(np.asarray(x, dtype=np.float64))
    xm = float(x_max if x_max is not None else x_arr.max())
    pos, salto = _saltos_psi(xm)
    if pos.size == 0:
        out = np.zeros_like(x_arr)
    else:
        cum = np.cumsum(salto)
        idx = np.searchsorted(pos, x_arr, side="right")
        out = np.where(idx > 0, cum[np.clip(idx - 1, 0, len(cum) - 1)], 0.0)
    return out if np.ndim(x) > 0 else float(out[0])


# ---------------------------------------------------------------------------
# Aproximación por fórmula explícita
# ---------------------------------------------------------------------------
def psi_aprox(x, gammas):
    """ψ_N(x) por la fórmula explícita truncada a los ceros `gammas` (partes
    imaginarias γ_n). Acepta x escalar o array."""
    x_arr = np.atleast_1d(np.asarray(x, dtype=np.float64))
    g = np.asarray(gammas, dtype=np.float64)[:, None]  # (N, 1)
    lnx = np.log(x_arr)[None, :]                        # (1, M)
    suma = np.real(np.exp(1j * g * lnx) / (0.5 + 1j * g)).sum(axis=0)  # (M,)
    out = (x_arr - 2.0 * np.sqrt(x_arr) * suma
           - LOG_2PI - 0.5 * np.log1p(-x_arr ** -2))
    return out if np.ndim(x) > 0 else float(out[0])


def error_medio(x_grid, gammas, psi_ex=None) -> float:
    """Error absoluto medio |ψ_N − ψ| sobre una malla (mide la convergencia con N)."""
    x_grid = np.asarray(x_grid, dtype=np.float64)
    if psi_ex is None:
        psi_ex = psi_exacta(x_grid, x_grid.max())
    return float(np.mean(np.abs(psi_aprox(x_grid, gammas) - psi_ex)))
