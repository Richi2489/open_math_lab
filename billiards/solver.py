"""
solver.py — espectro del Laplaciano (Helmholtz Dirichlet) en dominios 2D
========================================================================

Resuelve −∇²ψ = E ψ con ψ=0 en el borde (Dirichlet) por DIFERENCIAS FINITAS:
Laplaciano de 5 puntos sobre un grid regular, con una máscara que define el dominio
(los nodos fuera de la máscara imponen ψ=0). Los autovalores más bajos se obtienen con
`scipy.sparse.linalg.eigsh` en modo shift-invert (sigma=0).

Sin dependencias de mallado externo (no scikit-fem/gmsh): reproducible y cacheable. La
frontera "escalonada" del grid es una aproximación; por eso los modos bajos se desvían y
el match estadístico mejora ALTO en el espectro (se descartan los más bajos en el análisis).

Dominios incluidos:
  - rectángulo (también analítico, ver `eigs_rectangulo`);
  - cuarto de estadio de Bunimovich (desimetrizado), Dirichlet en TODO el borde del cuarto
    (incluidos los ejes de simetría) = un sector de paridad bien definido.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh


# ---------------------------------------------------------------------------
# Predicados de dominio (True dentro)
# ---------------------------------------------------------------------------
def rectangulo_dentro(X, Y, a, b):
    return (X > 0) & (X < a) & (Y > 0) & (Y < b)


def cuarto_estadio_dentro(X, Y, a=1.0, r=1.0):
    """Cuarto de estadio de Bunimovich: recta [0,a]×[0,r] + cuarto de disco radio r en x>a."""
    rect = (X > 0) & (X <= a) & (Y > 0) & (Y < r)
    arco = (X > a) & (Y > 0) & ((X - a) ** 2 + Y ** 2 < r ** 2)
    return rect | arco


# ---------------------------------------------------------------------------
# Solver de diferencias finitas
# ---------------------------------------------------------------------------
def resolver_dirichlet(dentro, x0, x1, y0, y1, h, n_modos):
    """Autovalores E (=k²) más bajos del Laplaciano Dirichlet en el dominio `dentro`.

    `dentro(X, Y)` → máscara booleana. Devuelve un array ordenado de `n_modos` autovalores.
    """
    xs = np.arange(x0 + h, x1, h)
    ys = np.arange(y0 + h, y1, h)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    mask = dentro(X, Y)
    nx, ny = X.shape

    ids = -np.ones((nx, ny), dtype=np.int64)
    ids[mask] = np.arange(int(mask.sum()))
    N = int(mask.sum())
    if N < n_modos + 5:
        raise ValueError(f"grid demasiado grueso: {N} nodos para {n_modos} modos.")

    inv_h2 = 1.0 / h ** 2
    filas = [ids[mask]]
    cols = [ids[mask]]
    datos = [np.full(N, 4.0 * inv_h2)]

    # vecinos en +x y +y (la simetría agrega los −x, −y)
    for sl_a, sl_b in (((slice(0, -1), slice(None)), (slice(1, None), slice(None))),
                       ((slice(None), slice(0, -1)), (slice(None), slice(1, None)))):
        m = mask[sl_a] & mask[sl_b]
        ia = ids[sl_a][m]
        ib = ids[sl_b][m]
        off = np.full(ia.size, -inv_h2)
        filas += [ia, ib]
        cols += [ib, ia]
        datos += [off, off]

    A = sp.coo_matrix((np.concatenate(datos),
                       (np.concatenate(filas), np.concatenate(cols))),
                      shape=(N, N)).tocsc()
    vals = eigsh(A, k=n_modos, sigma=0.0, which="LM", return_eigenvectors=False)
    return np.sort(vals.real)


# ---------------------------------------------------------------------------
# Rectángulo analítico (dominio integrable de control)
# ---------------------------------------------------------------------------
def eigs_rectangulo(a, b, n_modos):
    """Autovalores Dirichlet exactos del rectángulo [0,a]×[0,b]:
    E_{nm} = π²(n²/a² + m²/b²), n,m ≥ 1. Devuelve los `n_modos` más bajos, ordenados."""
    # cota de energía por Weyl, con margen
    E_max = (n_modos * 4 * np.pi / (a * b)) * 2.0 + 50.0
    n_top = int(np.sqrt(E_max) * a / np.pi) + 2
    m_top = int(np.sqrt(E_max) * b / np.pi) + 2
    n = np.arange(1, n_top + 1)[:, None]
    m = np.arange(1, m_top + 1)[None, :]
    E = np.pi ** 2 * ((n / a) ** 2 + (m / b) ** 2)
    E = np.sort(E.ravel())
    return E[:n_modos]
