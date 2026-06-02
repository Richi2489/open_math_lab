"""
engine.py — el mapa de Collatz
==============================

Define el mapa clásico de Collatz y su versión "acelerada" sobre los impares,
que es la que usa la heurística probabilística.

- Mapa clásico:  n par -> n/2 ; n impar -> 3n+1.
- Mapa acelerado (sobre impares):  n -> (3n+1) / 2^v,  donde v = v_2(3n+1).
  Como 3n+1 es par cuando n es impar, siempre v >= 1, y el resultado vuelve a ser impar.

Las funciones escalares usan enteros de Python (precisión arbitraria) para que las
trayectorias largas NO desborden. Las funciones vectorizadas usan numpy int64 y solo
deben llamarse con valores acotados (ver `paso_acelerado_vector`).

`numba` es opcional: si está instalado se usa para acelerar el conteo de tiempos de
parada sobre valores que caben en int64.
"""

from __future__ import annotations

import numpy as np

# --- numba opcional --------------------------------------------------------
try:
    from numba import njit

    _TIENE_NUMBA = True
except ImportError:  # pragma: no cover - depende del entorno
    _TIENE_NUMBA = False


# ---------------------------------------------------------------------------
# Mapa clásico
# ---------------------------------------------------------------------------
def trayectoria(n: int) -> list[int]:
    """Trayectoria completa de `n` bajo el mapa clásico hasta llegar a 1.

    Devuelve la lista [n, ..., 1] (incluye el punto inicial y el 1 final).
    """
    if n < 1:
        raise ValueError("n debe ser un entero positivo.")
    seq = [n]
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        seq.append(n)
    return seq


def tiempo_total_parada(n: int) -> int:
    """Tiempo total de parada: número de pasos del mapa clásico para llegar a 1."""
    if n < 1:
        raise ValueError("n debe ser un entero positivo.")
    pasos = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else 3 * n + 1
        pasos += 1
    return pasos


# ---------------------------------------------------------------------------
# Valuación 2-ádica
# ---------------------------------------------------------------------------
def valuacion_2(m: int) -> int:
    """v_2(m): número de veces que 2 divide a m (orden 2-ádico).

    v_2(0) se define aquí como 0 por conveniencia (no debería ocurrir en Collatz).
    """
    if m == 0:
        return 0
    v = 0
    while m & 1 == 0:
        m >>= 1
        v += 1
    return v


# ---------------------------------------------------------------------------
# Mapa acelerado sobre impares
# ---------------------------------------------------------------------------
def paso_acelerado(n: int) -> tuple[int, int]:
    """Un paso del mapa acelerado sobre un impar `n`.

    Calcula m = 3n+1, v = v_2(m), y devuelve (m / 2^v, v).
    El primer elemento del par siempre es impar.
    """
    m = 3 * n + 1
    v = valuacion_2(m)
    return m >> v, v


def secuencia_v(n: int) -> np.ndarray:
    """Secuencia de v's a lo largo de la trayectoria acelerada de `n` hasta llegar a 1.

    Si `n` es par, primero se reduce a impar (se descartan esos factores de 2,
    que no aportan información sobre la dinámica del mapa acelerado).
    """
    n = int(n)
    if n < 1:
        raise ValueError("n debe ser un entero positivo.")
    if n % 2 == 0:
        n >>= valuacion_2(n)
    vs: list[int] = []
    while n != 1:
        n, v = paso_acelerado(n)
        vs.append(v)
    return np.asarray(vs, dtype=np.int64)


# ---------------------------------------------------------------------------
# Versiones vectorizadas (numpy int64). Usar solo con valores acotados.
# ---------------------------------------------------------------------------
def valuacion_2_vector(m: np.ndarray) -> np.ndarray:
    """v_2 vectorizado para un array de enteros int64 (todos > 0)."""
    m = np.asarray(m, dtype=np.int64).copy()
    v = np.zeros_like(m)
    activos = (m & 1) == 0
    while activos.any():
        m[activos] >>= 1
        v[activos] += 1
        activos = (m & 1) == 0
    return v


def paso_acelerado_vector(n: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Un paso del mapa acelerado para un array de impares (int64).

    Devuelve (siguiente, v). Pensado para muestreo estadístico: con `n` impar
    hasta ~2^40, el cálculo 3n+1 cabe holgadamente en int64 sin desbordar.
    """
    n = np.asarray(n, dtype=np.int64)
    m = 3 * n + 1
    v = valuacion_2_vector(m)
    siguiente = m >> v
    return siguiente, v


# ---------------------------------------------------------------------------
# Aceleración opcional con numba (solo valores que caben en int64)
# ---------------------------------------------------------------------------
if _TIENE_NUMBA:  # pragma: no cover - depende del entorno

    @njit(cache=True)
    def _tiempo_total_parada_int64(n: np.int64) -> np.int64:
        pasos = 0
        while n != 1:
            if n % 2 == 0:
                n = n // 2
            else:
                n = 3 * n + 1
            pasos += 1
        return pasos

    def tiempo_total_parada_rapido(n: int) -> int:
        """Versión acelerada con numba. CAVEAT: usa int64, válida solo si toda la
        trayectoria de `n` cabe en int64 (no hay protección contra desbordamiento)."""
        return int(_tiempo_total_parada_int64(np.int64(n)))

else:

    def tiempo_total_parada_rapido(n: int) -> int:
        """numba no está instalado: cae a la implementación pura de Python."""
        return tiempo_total_parada(n)
