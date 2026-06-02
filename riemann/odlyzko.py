"""
odlyzko.py — ingesta de tablas de ceros de Odlyzko (alta altura T)
==================================================================

Carga ceros de ζ desde un archivo LOCAL de las tablas de Odlyzko (no descarga en
runtime). El dataset objetivo de la iteración 2 es `zeros6`: los primeros 2,001,052
ceros, precisión ~4×10⁻⁹, altura hasta T ≈ 1.13×10⁶. Con ~2M de espaciados el piso de
ruido del KS colapsa de ~0.06 (con 500) a ~0.001, así que la corrección conocida de
altura finita al GUE se vuelve visible — NO es una anomalía.

El loader devuelve el MISMO tipo de array (np.float64 1D de partes imaginarias γ) que
`riemann.zeros`, para que el resto del lab no distinga la fuente.

Formatos aceptados:
  - txt de una columna (con o sin espacios al inicio);
  - csv de una columna;
  - csv con cabecera y columnas nombradas (se detecta gamma/zero/zeros/t/ordinate).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

_NOMBRES_COL = ("gamma", "zero", "zeros", "t", "ordinate", "im", "imag")


def _es_cabecera(linea: str) -> bool:
    """True si la primera celda de la línea NO es un número (=> es cabecera)."""
    primera = linea.replace(",", " ").split()[0]
    try:
        float(primera)
        return False
    except ValueError:
        return True


def cargar_odlyzko(path) -> np.ndarray:
    """Carga las partes imaginarias γ de un archivo de ceros de Odlyzko.

    Devuelve un np.ndarray 1D float64 (mismo formato que `riemann.zeros.cargar_ceros`).
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"No existe {path}. Obtén el dataset de Odlyzko (p.ej. zeros6) y colócalo ahí; "
            "el script NO descarga en runtime a propósito (evita fragilidad)."
        )
    with open(path, "r", encoding="utf-8", errors="strict") as fh:
        primera = ""
        for linea in fh:
            if linea.strip():
                primera = linea
                break

    if not primera:
        raise ValueError(f"{path} está vacío.")

    coma = "," in primera
    col = 0
    saltar = 0
    if _es_cabecera(primera):
        saltar = 1
        sep = "," if coma else None
        cabeceras = [c.strip().lower() for c in
                     (primera.split(",") if coma else primera.split())]
        col = next((i for i, c in enumerate(cabeceras)
                    if any(n == c or n in c for n in _NOMBRES_COL)), 0)

    delim = "," if coma else None
    datos = np.loadtxt(path, delimiter=delim, skiprows=saltar, usecols=col, dtype=np.float64)
    return np.atleast_1d(datos).astype(np.float64)


def validar_zeros(gammas: np.ndarray, min_len: int = 2) -> None:
    """Valida el array de ceros. Lanza ValueError ante cualquier problema."""
    g = np.asarray(gammas)
    if g.ndim != 1:
        raise ValueError("se esperaba un array 1D.")
    if not np.issubdtype(g.dtype, np.number):
        raise ValueError("el array no es numérico.")
    if g.size < min_len:
        raise ValueError(f"muy pocos ceros: {g.size} < {min_len}.")
    if not np.all(np.isfinite(g)):
        raise ValueError("hay NaN o inf.")
    d = np.diff(g)
    if np.any(d == 0):
        raise ValueError("hay ceros duplicados (spacing nulo).")
    if np.any(d < 0):
        raise ValueError("los ceros no son estrictamente ascendentes (spacing negativo).")


def resumir_dataset(gammas: np.ndarray) -> dict:
    """Resumen del dataset: n, rango de altura T y spacing crudo medio."""
    g = np.asarray(gammas, dtype=np.float64)
    return {
        "n": int(g.size),
        "T_min": float(g[0]),
        "T_max": float(g[-1]),
        "spacing_crudo_medio": float(np.diff(g).mean()),
    }


def tomar_ventana(gammas, start: int = 0, stop=None, n=None) -> np.ndarray:
    """Subconjunto contiguo de ceros: [start:stop] o `n` a partir de `start`."""
    g = np.asarray(gammas, dtype=np.float64)
    if n is not None:
        return g[start:start + n].copy()
    return g[start:stop].copy()


def dividir_en_ventanas(gammas, k: int) -> list:
    """Parte los ceros en `k` ventanas contiguas NO solapadas (por altura)."""
    if k < 1:
        raise ValueError("k debe ser >= 1.")
    return [w.copy() for w in np.array_split(np.asarray(gammas, dtype=np.float64), k)]
