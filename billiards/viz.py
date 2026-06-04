"""
viz.py — figuras del laboratorio de billares
=============================================

Backend 'Agg'. Reusa la maquinaria de histogramas de `riemann.metrics` y las curvas
teóricas de `riemann.gue`.
"""

from __future__ import annotations

import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from riemann import gue, metrics  # noqa: E402


def _asegurar_dir(ruta):
    carpeta = os.path.dirname(ruta)
    if carpeta:
        os.makedirs(carpeta, exist_ok=True)


def par_espaciados(s_rect, s_stad, ruta) -> str:
    """LA figura clave: las dos distribuciones de espaciado lado a lado.
    Izquierda: rectángulo (integrable → Poisson, sin hoyo en s→0).
    Derecha: estadio (caótico → GOE, con hoyo de repulsión). Mismo solver, sólo cambió
    la geometría de la mesa."""
    _asegurar_dir(ruta)
    xx = np.linspace(0, 4, 300)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    c, d = metrics.histograma(s_rect)
    a1.bar(c, d, width=(c[1] - c[0]) * 0.9, color="#c0a06a", alpha=0.7, label="rectángulo")
    a1.plot(xx, np.exp(-xx), "-", color="#a5483b", lw=2, label="Poisson  e⁻ˢ")
    a1.set_title("Integrable (rectángulo) → Poisson\nsin repulsión: P(s) máximo en s→0")
    a1.set_xlabel("espaciado normalizado s"); a1.set_ylabel("P(s)")
    a1.legend(); a1.grid(True, alpha=0.3); a1.set_xlim(0, 4)

    c, d = metrics.histograma(s_stad)
    a2.bar(c, d, width=(c[1] - c[0]) * 0.9, color="#5a86b4", alpha=0.7, label="estadio")
    a2.plot(xx, gue.wigner_goe(xx), "-", color="#23496e", lw=2, label="GOE (Wigner)")
    a2.plot(xx, np.exp(-xx), ":", color="#a5483b", lw=1.2, label="Poisson (ref.)")
    a2.set_title("Caótico (estadio) → GOE\nhoyo de repulsión: P(s)→0 cuando s→0")
    a2.set_xlabel("espaciado normalizado s")
    a2.legend(); a2.grid(True, alpha=0.3); a2.set_xlim(0, 4)

    fig.suptitle("Caos vs integrabilidad: mismo solver, sólo cambia la mesa (Bohigas–Giannoni–Schmit)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def cdf_par(s_rect, s_stad, s_goe, s_pois, ruta) -> str:
    """CDF empíricas: rectángulo y estadio frente a Poisson y GOE."""
    _asegurar_dir(ruta)
    fig, ax = plt.subplots(figsize=(9, 5))
    for s, et, col, ls in [(s_pois, "Poisson", "#a5483b", ":"),
                           (s_goe, "GOE", "#23496e", ":"),
                           (s_rect, "rectángulo", "#c0a06a", "-"),
                           (s_stad, "estadio", "#5a86b4", "-")]:
        x, F = metrics.cdf_empirica(s)
        ax.plot(x, F, ls, color=col, lw=1.6, label=et)
    ax.set_xlim(0, 4); ax.set_xlabel("espaciado normalizado s"); ax.set_ylabel("CDF")
    ax.set_title("CDF de espaciados: el rectángulo sigue a Poisson, el estadio a GOE")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def conteo_weyl(E, area, perim, ruta, etiqueta="dominio") -> str:
    """Escalera real N(E) vs el conteo suave de Weyl (área + perímetro)."""
    from . import weyl
    _asegurar_dir(ruta)
    E = np.sort(np.asarray(E, dtype=np.float64))
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.step(E, np.arange(1, len(E) + 1), where="post", color="#5a86b4", lw=1.1,
            label=f"escalera N(E) ({etiqueta})")
    EE = np.linspace(E[0], E[-1], 400)
    ax.plot(EE, weyl.N_weyl(EE, area, perim), "-", color="#a5483b", lw=1.8,
            label="Weyl (A/4π)E − (P/4π)√E")
    ax.set_xlabel("E"); ax.set_ylabel("N(E)")
    ax.set_title(f"Ley de Weyl: {etiqueta}")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta
