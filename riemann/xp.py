"""
xp.py — el operador xp de Berry–Keating: tocar la silueta, no el rostro
=======================================================================

Demostración controlada (NO es un intento de RH) de por qué el Hamiltoniano `H = xp`
de Berry–Keating "casi funciona" y por qué no basta.

Dos hechos:

1. SILUETA (densidad suave). El conteo semiclásico de fase del Hamiltoniano CLÁSICO
   H = xp, con la regularización de Berry–Keating, es

       N(E) = (E/2π)·(log(E/2π) − 1) + 7/8 + …

   que es EXACTAMENTE la parte suave de Riemann–von Mangoldt ⟨N(T)⟩ — la misma densidad
   que usa el unfolding de este repo (`spacing.N_suave`). Por eso xp es "el sospechoso".

2. ROSTRO (espectro fino). El generador de dilatación auto-adjunto
   H = −i(x d/dx + 1/2) es, en la variable u = log x (con la medida unitaria correcta),
   el operador de momento −i d/du. En una caja u ∈ [0, L] su espectro es un PEINE casi
   uniforme E_n ≈ 2πn/L: rígido, sin las correcciones que distinguen a los ceros t_n, y
   sin repulsión tipo GUE. Captura la densidad media, no la estructura fina.
"""

from __future__ import annotations

import numpy as np

from . import spacing

DOS_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# 1. Silueta: conteo semiclásico de xp
# ---------------------------------------------------------------------------
def conteo_semiclasico_xp(E):
    """N(E) semiclásico (Berry–Keating, regularizado) de H = xp:

        N(E) = (E/2π)(log(E/2π) − 1) + 7/8.

    Es IDÉNTICO a la parte suave de Riemann–von Mangoldt (`spacing.N_suave`): esa es la
    coincidencia que motiva todo el programa espectral.
    """
    return spacing.N_suave(np.asarray(E, dtype=np.float64))


def conteo_leading_order(E):
    """Forma a leading order (E/2π)(log E − 1), que coincide con la anterior salvo los
    términos sub-dominantes (el −log 2π y el 7/8)."""
    E = np.asarray(E, dtype=np.float64)
    return (E / DOS_PI) * (np.log(E) - 1.0)


# ---------------------------------------------------------------------------
# 2. Rostro: el generador de dilatación discretizado
# ---------------------------------------------------------------------------
def operador_dilatacion(N: int, L: float) -> np.ndarray:
    """Discretiza H = −i d/du (el generador de dilatación en u = log x) en una caja
    periódica u ∈ [0, L), por método espectral de Fourier. Devuelve la matriz Hermitiana
    N×N (N par). Sus autovalores forman el peine 2π·m/L.

    El término +1/2 de H = −i(x d/dx + 1/2) queda absorbido por la transformación unitaria
    a la variable u (la medida dx = e^u du), así que en u el operador es exactamente el
    momento −i d/du.
    """
    if N % 2 != 0:
        raise ValueError("N debe ser par para la matriz espectral de Fourier.")
    h = DOS_PI / N
    D = np.zeros((N, N), dtype=np.float64)  # matriz de derivación espectral (periodo 2π)
    for i in range(N):
        for j in range(N):
            if i != j:
                d = i - j
                D[i, j] = 0.5 * ((-1.0) ** d) / np.tan(d * h / 2.0)
    D *= DOS_PI / L  # escalar al periodo L: d/du
    return -1j * D  # Hermitiana


def espectro_xp(N: int, L: float) -> np.ndarray:
    """Autovalores (ordenados) del generador de dilatación en la caja [0, L)."""
    return np.sort(np.linalg.eigvalsh(operador_dilatacion(N, L)))


def niveles_positivos(E: np.ndarray, recorte: float = 0.5) -> np.ndarray:
    """Autovalores positivos en el bulk (descarta una fracción `recorte` del extremo
    superior, donde la rejilla introduce dispersión). El peine resultante es ~uniforme."""
    E = np.sort(np.asarray(E, dtype=np.float64))
    pos = E[E > 0]
    if recorte > 0 and len(pos) > 4:
        pos = pos[: int(len(pos) * (1.0 - recorte))]
    return pos


def espaciados_peine(niveles: np.ndarray) -> np.ndarray:
    """Espaciados del peine xp, unfolded por su densidad media (constante): s_n con
    media 1. Para un peine perfectamente uniforme, todos los s_n ≈ 1 (rígido)."""
    gaps = np.diff(np.sort(np.asarray(niveles, dtype=np.float64)))
    return gaps / gaps.mean()
