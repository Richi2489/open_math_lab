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
import numpy as np

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


# ===========================================================================
# Figuras de la iteración 2
# ===========================================================================
def acf_con_banda(resultado, clave, ruta, titulo) -> str:
    """ACF (bruto o corregido) vs lag, con banda analítica ±1.96/sqrt(N_eff).

    `resultado`: dict de acf.acf_agregado. `clave`: 'rho_bruto' o 'rho_corregido'.
    """
    _asegurar_dir(ruta)
    lags = resultado["lags"][1:]
    rho = resultado[clave][1:]
    banda = resultado["banda"][1:]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.fill_between(lags, -banda, banda, color="#cccccc", alpha=0.6,
                    label="banda ±1.96/√N_eff")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.plot(lags, rho, "o-", color="#3b6ea5", lw=1.5, label=clave.replace("_", " "))
    ax.set_xlabel("lag")
    ax.set_ylabel("autocorrelación")
    ax.set_title(titulo)
    ax.set_xticks(lags)
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def lag1_por_grupo(etiquetas, valores, err_inf, err_sup, ruta, titulo, xlabel) -> str:
    """Barras del lag-1 corregido por grupo (bucket de L o banda de bits) con barras
    de error asimétricas (p.ej. percentiles de permutación)."""
    _asegurar_dir(ruta)
    x = range(len(etiquetas))
    valores = np.asarray(valores)
    yerr = np.vstack([np.asarray(err_inf), np.asarray(err_sup)])
    fig, ax = plt.subplots(figsize=(9, 5))
    colores = ["#3b6ea5" if v > 0 else "#a5483b" for v in valores]
    ax.bar(list(x), valores, color=colores, alpha=0.85)
    ax.errorbar(list(x), valores, yerr=yerr, fmt="none", ecolor="black", capsize=4, lw=1)
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(etiquetas, rotation=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("lag-1 corregido")
    ax.set_title(titulo)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def histograma_longitudes(longitudes, ruta="outputs/it2_longitudes.png") -> str:
    """Histograma de los largos L de las trayectorias muestreadas."""
    _asegurar_dir(ruta)
    longitudes = np.asarray(longitudes)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.hist(longitudes, bins=50, color="#5a9367", alpha=0.85, edgecolor="white", lw=0.3)
    ax.axvline(float(longitudes.mean()), color="crimson", ls="--", lw=1.5,
               label=f"media = {longitudes.mean():.1f}")
    ax.set_xlabel("L (número de pasos acelerados)")
    ax.set_ylabel("frecuencia")
    ax.set_title("Distribución de largos de trayectoria")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def distribucion_v(v_concat, ruta="outputs/it2_dist_v.png", k_max=12) -> str:
    """Distribución empírica de v frente al modelo geométrico P(v=k)=2^-k."""
    _asegurar_dir(ruta)
    v_concat = np.asarray(v_concat)
    ks = np.arange(1, k_max + 1)
    emp = np.array([(v_concat == k).mean() for k in ks])
    teo = 2.0 ** (-ks.astype(np.float64))
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(ks - 0.18, emp, width=0.36, color="#3b6ea5", alpha=0.85, label="empírico")
    ax.bar(ks + 0.18, teo, width=0.36, color="#d6a04c", alpha=0.85, label="2^-k (modelo)")
    ax.set_yscale("log")
    ax.set_xlabel("v")
    ax.set_ylabel("P(v) (escala log)")
    ax.set_title("Distribución marginal de v vs modelo geométrico")
    ax.set_xticks(ks)
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta
