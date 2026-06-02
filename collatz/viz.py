"""
viz.py — gráficas del laboratorio
=================================

Tres figuras:
  - trayectorias en escala logarítmica,
  - histograma de tiempos de parada,
  - curva de autocorrelación de los v's.

Usa el backend 'Agg' (sin ventana) para poder guardar archivos en entornos sin GUI,
como una terminal o CI. Todas las funciones devuelven la ruta del archivo guardado.
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")  # backend sin ventana: solo guarda a disco
import matplotlib.pyplot as plt  # noqa: E402

from .engine import trayectoria  # noqa: E402


def _asegurar_dir(ruta: str) -> None:
    carpeta = os.path.dirname(ruta)
    if carpeta:
        os.makedirs(carpeta, exist_ok=True)


def grafica_trayectorias(semillas, ruta="outputs/trayectorias.png") -> str:
    """Dibuja las trayectorias (escala log en el eje Y) de las semillas dadas."""
    _asegurar_dir(ruta)
    fig, ax = plt.subplots(figsize=(9, 5))
    for n in semillas:
        seq = trayectoria(int(n))
        ax.plot(range(len(seq)), seq, lw=1, alpha=0.8, label=f"n={n}")
    ax.set_yscale("log")
    ax.set_xlabel("paso")
    ax.set_ylabel("valor (escala log)")
    ax.set_title("Trayectorias de Collatz")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def histograma_tiempos_parada(tiempos, ruta="outputs/tiempos_parada.png") -> str:
    """Histograma de la distribución de tiempos totales de parada."""
    _asegurar_dir(ruta)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(tiempos, bins=60, color="#3b6ea5", alpha=0.85, edgecolor="white", lw=0.3)
    ax.axvline(
        float(tiempos.mean()),
        color="crimson",
        ls="--",
        lw=1.5,
        label=f"media = {tiempos.mean():.1f}",
    )
    ax.set_xlabel("tiempo total de parada")
    ax.set_ylabel("frecuencia")
    ax.set_title("Distribución de tiempos de parada")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def grafica_autocorrelacion(lags, acf, ruta="outputs/autocorrelacion.png") -> str:
    """Diagrama de tallo de la autocorrelación de los v's frente al lag."""
    _asegurar_dir(ruta)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.stem(lags, acf)
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xlabel("lag")
    ax.set_ylabel("autocorrelación")
    ax.set_title("Autocorrelación de los v's (independencia asumida = 0 en lag>=1)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta
