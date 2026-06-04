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


def cdf_diferencia(s_zeta, s_gue, ruta, x_max=4.0) -> str:
    """Diferencia de CDF empíricas F_ζ(x) − F_GUE(x) (la 'firma' del residuo)."""
    _asegurar_dir(ruta)
    x = np.linspace(0, x_max, 400)
    sz = np.sort(np.asarray(s_zeta, dtype=np.float64))
    sg = np.sort(np.asarray(s_gue, dtype=np.float64))
    Fz = np.searchsorted(sz, x, side="right") / len(sz)
    Fg = np.searchsorted(sg, x, side="right") / len(sg)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x, Fz - Fg, color="#3b6ea5", lw=1.5)
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xlabel("espaciado normalizado s")
    ax.set_ylabel("F_ζ(s) − F_GUE(s)")
    ax.set_title("Diferencia de CDF (ζ Odlyzko − GUE)")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def metricas_por_ventana(centros_T, ks, wass, ruta) -> str:
    """KS y Wasserstein (ζ vs GUE) por ventana de altura."""
    _asegurar_dir(ruta)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.5))
    a1.plot(centros_T, ks, "o-", color="#3b6ea5", lw=1.5)
    a1.set_xlabel("altura T (centro de ventana)")
    a1.set_ylabel("KS(ζ vs GUE)")
    a1.set_title("KS por ventana de altura")
    a1.grid(True, alpha=0.3)
    a2.plot(centros_T, wass, "o-", color="#5a9367", lw=1.5)
    a2.set_xlabel("altura T (centro de ventana)")
    a2.set_ylabel("Wasserstein(ζ vs GUE)")
    a2.set_title("Wasserstein por ventana de altura")
    a2.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def comparar_datasets(etiquetas, valores, ruta, ylabel="Wasserstein(ζ vs GUE)") -> str:
    """Barras de una métrica de desviación a GUE a través de datasets/ventanas
    (p.ej. baja altura mpmath vs Odlyzko alta-T)."""
    _asegurar_dir(ruta)
    x = range(len(etiquetas))
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(list(x), valores, color="#3b6ea5", alpha=0.85)
    ax.set_xticks(list(x))
    ax.set_xticklabels(etiquetas, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title("Desviación a GUE: baja altura vs Odlyzko alta-T")
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def formula_explicita_panel(x, psi_ex, curvas_por_N, ruta) -> str:
    """La figura de la §4: ψ(x) exacta (escalera de primos) vs su reconstrucción por la
    fórmula explícita con N ceros crecientes. Panel 2×2."""
    _asegurar_dir(ruta)
    x = np.asarray(x, dtype=np.float64)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    for ax, (N, ap) in zip(axes.ravel(), curvas_por_N):
        ax.plot(x, psi_ex, color="black", lw=1.3, label="ψ(x) exacta (primos)")
        ax.plot(x, ap, color="#3b6ea5", lw=1.0, label=f"fórmula explícita, N={N}")
        ax.set_title(f"N = {N} ceros")
        ax.set_ylabel("ψ(x)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc="upper left")
    for ax in axes[-1]:
        ax.set_xlabel("x")
    fig.suptitle("Reconstrucción del conteo de primos ψ(x) desde los ceros de ζ "
                 "(fórmula explícita de Riemann)", fontsize=12)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def formula_explicita_gibbs(x, psi_ex, curvas_por_N, ruta) -> str:
    """Zoom cerca de unos saltos para mostrar el fenómeno de Gibbs (overshoot) como
    comportamiento esperado de la suma truncada, no como error."""
    _asegurar_dir(ruta)
    x = np.asarray(x, dtype=np.float64)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x, psi_ex, color="black", lw=1.6, label="ψ(x) exacta")
    colores = ["#ccdcea", "#9bbcd8", "#5a86b4", "#23496e"]
    for (N, ap), c in zip(curvas_por_N, colores):
        ax.plot(x, ap, color=c, lw=1.1, label=f"N={N}")
    ax.set_xlabel("x")
    ax.set_ylabel("ψ(x)")
    ax.set_title("Fenómeno de Gibbs cerca de los saltos (esperado, no es error)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(ruta, dpi=130)
    plt.close(fig)
    return ruta


def silueta_conteo(t_zeros, ruta, n_zeros=400) -> str:
    """A.1: escalera real N(t)=#{t_n≤t} de los ceros vs el conteo semiclásico de xp
    (= parte suave de Riemann–von Mangoldt). Coinciden en promedio → la SILUETA."""
    from . import spacing
    _asegurar_dir(ruta)
    t = np.sort(np.asarray(t_zeros, dtype=np.float64))[:n_zeros]
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.step(t, np.arange(1, len(t) + 1), where="post", color="#3b6ea5", lw=1.2,
            label="escalera real N(t) de los ceros")
    tt = np.linspace(t[0], t[-1], 400)
    ax.plot(tt, spacing.N_suave(tt), "-", color="#a5483b", lw=1.8,
            label="conteo semiclásico de xp = N_suave(T)")
    ax.set_xlabel("T"); ax.set_ylabel("N(T)")
    ax.set_title("A.1 — xp reproduce la densidad suave (silueta)")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(ruta, dpi=130); plt.close(fig)
    return ruta


def peine_xp(niveles, L, ruta) -> str:
    """A.2: los autovalores positivos del generador de dilatación discretizado:
    un peine ~uniforme E_n ≈ 2πn/L."""
    _asegurar_dir(ruta)
    niveles = np.sort(np.asarray(niveles, dtype=np.float64))
    n = np.arange(1, len(niveles) + 1)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(n, niveles, "o", color="#5a9367", ms=3, label="autovalores xp (caja)")
    ax.plot(n, 2 * np.pi * n / L, "-", color="#999999", lw=1,
            label="2πn/L (peine uniforme)")
    ax.set_xlabel("n"); ax.set_ylabel("E_n")
    ax.set_title("A.2 — espectro del generador de dilatación: peine uniforme")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(ruta, dpi=130); plt.close(fig)
    return ruta


def niveles_vs_zeros(t_zeros, ruta, n_mostrar=400) -> str:
    """A.3a: t_n reales vs un peine uniforme escalado al mismo rango. El peine es la
    cuerda recta; los ceros se curvan (densidad creciente) → DIVERGEN."""
    _asegurar_dir(ruta)
    t = np.sort(np.asarray(t_zeros, dtype=np.float64))[:n_mostrar]
    n = np.arange(1, len(t) + 1)
    peine = t[0] + (t[-1] - t[0]) * (n - 1) / (len(t) - 1)  # comb lineal mismo rango
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(n, t, "-", color="#3b6ea5", lw=1.5, label="ceros de ζ: t_n")
    ax.plot(n, peine, "--", color="#5a9367", lw=1.5, label="peine xp uniforme (mismo rango)")
    ax.set_xlabel("n"); ax.set_ylabel("valor")
    ax.set_title("A.3 — niveles individuales: xp (recto) vs t_n (curvado) divergen")
    ax.legend(); ax.grid(True, alpha=0.3)
    fig.tight_layout(); fig.savefig(ruta, dpi=130); plt.close(fig)
    return ruta


def espaciados_xp(s_xp, s_gue, s_zeta, ruta) -> str:
    """A.3b: distribución de espaciados. xp = pico rígido en s≈1 (picket-fence); GUE y
    ζ siguen la sorpresa de Wigner con repulsión. xp NO genera estadística GUE."""
    from . import gue, metrics
    _asegurar_dir(ruta)
    fig, ax = plt.subplots(figsize=(9, 5))
    cg, dg = metrics.histograma(s_gue)
    ax.plot(cg, dg, "-", color="#5a9367", lw=1.4, label="GUE")
    cz, dz = metrics.histograma(s_zeta)
    ax.plot(cz, dz, "o-", color="#3b6ea5", ms=3, lw=1.1, label="ζ (ceros)")
    xx = np.linspace(0, 4, 300)
    ax.plot(xx, gue.wigner_gue(xx), ":", color="#a5483b", lw=1.5, label="Wigner GUE")
    ax.axvline(float(np.mean(s_xp)), color="#d6a04c", lw=2.2,
               label="xp: pico rígido en s≈1 (sin repulsión)")
    ax.set_xlabel("espaciado normalizado s"); ax.set_ylabel("densidad P(s)")
    ax.set_title("A.3 — estadística de espaciados: xp ≠ GUE")
    ax.legend(); ax.grid(True, alpha=0.3); ax.set_xlim(0, 4)
    fig.tight_layout(); fig.savefig(ruta, dpi=130); plt.close(fig)
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
