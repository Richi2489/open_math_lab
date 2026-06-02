"""
viz.py — figuras del Riemann Statistical Lab
============================================

Backend 'Agg' (sin ventana) para guardar a disco en cualquier entorno.
"""

from __future__ import annotations

import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from . import gue, metrics  # noqa: E402


def _asegurar_dir(ruta):
    carpeta = os.path.dirname(ruta)
    if carpeta:
        os.makedirs(carpeta, exist_ok=True)


def espaciados_nn(s_zeta, s_gue, s_poisson, ruta="outputs/riemann_spacing.png") -> str:
    """Histograma de espaciados de vecino más cercano: ζ vs GUE vs Poisson + Wigner."""
    _asegurar_dir(ruta)
    fig, ax = plt.subplots(figsize=(9, 5))
    cz, dz = metrics.histograma(s_zeta)
    ax.plot(cz, dz, "o-", color="#3b6ea5", ms=3, lw=1.3, label="ζ (unfolded)")
    cg, dg = metrics.histograma(s_gue)
    ax.plot(cg, dg, "-", color="#5a9367", lw=1.3, label="GUE finito")
    cp, dp = metrics.histograma(s_poisson)
    ax.plot(cp, dp, "--", color="#999999", lw=1.2, label="Poisson")
    xx = np.linspace(0, 4, 300)
    ax.plot(xx, gue.wigner_gue(xx), ":", color="#a5483b", lw=1.6, label="Wigner GUE (teoría)")
    ax.set_xlabel("espaciado normalizado s")
    ax.set_ylabel("densidad P(s)")
    ax.set_title("Espaciados de vecino más cercano")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def cdf_espaciados(s_zeta, s_gue, s_poisson, ruta="outputs/riemann_cdf.png") -> str:
    """CDF empírica de los espaciados: ζ vs GUE vs Poisson."""
    _asegurar_dir(ruta)
    fig, ax = plt.subplots(figsize=(9, 5))
    for s, est, col in [(s_zeta, "ζ", "#3b6ea5"), (s_gue, "GUE", "#5a9367"),
                        (s_poisson, "Poisson", "#999999")]:
        x, F = metrics.cdf_empirica(s)
        ax.plot(x, F, lw=1.4, label=est, color=col)
    ax.set_xlim(0, 4)
    ax.set_xlabel("espaciado normalizado s")
    ax.set_ylabel("CDF")
    ax.set_title("CDF de espaciados")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def pares_r2(centros, r2_emp, ruta="outputs/riemann_r2.png") -> str:
    """Correlación de pares empírica de ζ vs la curva GUE 1 − (sin(πr)/(πr))²."""
    _asegurar_dir(ruta)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(centros, r2_emp, "o-", color="#3b6ea5", ms=3, lw=1.2, label="R₂ empírico (ζ)")
    rr = np.linspace(centros.min(), centros.max(), 400)
    ax.plot(rr, metrics.r2_teorico_gue(rr), "-", color="#a5483b", lw=1.6,
            label="1 − (sin(πr)/(πr))²  (GUE)")
    ax.axhline(1.0, color="#999999", ls="--", lw=1, label="Poisson (R₂=1)")
    ax.set_xlabel("r (en unidades de espaciado medio)")
    ax.set_ylabel("R₂(r)")
    ax.set_title("Correlación de pares (firma Montgomery–Odlyzko)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def tendencia_altura(centros_T, var_gap, ks, ruta="outputs/riemann_tendencia_altura.png") -> str:
    """Dos paneles: |brecha de varianza ζ−GUE| y KS(ζ vs GUE) por cuartil de altura.
    Si la brecha encoge al subir T, es efecto de altura finita."""
    _asegurar_dir(ruta)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5))
    a1.plot(centros_T, var_gap, "o-", color="#3b6ea5", lw=1.6)
    a1.set_xlabel("altura T (centro de cuartil)")
    a1.set_ylabel("|var(ζ) − var(GUE)|")
    a1.set_title("Brecha de varianza vs altura")
    a1.grid(True, alpha=0.3)
    a2.plot(centros_T, ks, "o-", color="#5a9367", lw=1.6)
    a2.set_xlabel("altura T (centro de cuartil)")
    a2.set_ylabel("KS(ζ vs GUE)")
    a2.set_title("KS vs altura")
    a2.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def estabilidad_altura(estab, ruta="outputs/riemann_estabilidad.png") -> str:
    """KS(ζ vs GUE) por ventana de altura T (¿estable o efecto de tamaño finito?)."""
    _asegurar_dir(ruta)
    T = [0.5 * (e["T_min"] + e["T_max"]) for e in estab]
    ks = [e["ks"] for e in estab]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(T, ks, "o-", color="#3b6ea5", lw=1.4)
    ax.set_xlabel("altura T (centro de ventana)")
    ax.set_ylabel("KS(ζ vs GUE finito)")
    ax.set_title("Estabilidad del match GUE por altura")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta
