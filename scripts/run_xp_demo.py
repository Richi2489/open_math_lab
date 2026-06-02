"""
run_xp_demo.py — "Touching the Wall": por qué el operador xp casi funciona
==========================================================================

Cierre conceptual del Riemann Statistical Lab. NO es un intento de RH. Es una
demostración controlada de la idea de Berry–Keating: el Hamiltoniano H = xp reproduce la
DENSIDAD MEDIA de los ceros (la silueta), pero no su espectro fino, ni la repulsión tipo
GUE, ni la huella aritmética de los primos (el rostro).

Escribe docs/riemann_missing_operator.md y figuras en docs/figures/riemann_xp/.
Lenguaje prudente.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, ValueError):  # pragma: no cover
    pass

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from riemann import gue, metrics, spacing, viz, xp, zeros  # noqa: E402

FIGDIR = "docs/figures/riemann_xp"
N_GRID = 600       # puntos del grid (matriz del generador de dilatación)
L_CAJA = 30.0      # longitud de la caja en u = log x
N_ZEROS = 2000     # ceros cacheados (mpmath) para comparar
SEED = 2026


def f(x, n=4):
    return f"{x:.{n}f}"


def main():
    R = []

    def linea(s=""):
        print(s)
        R.append(s)

    rng = np.random.default_rng(SEED)

    linea("=" * 72)
    linea("  Riemann Lab — \"Touching the Wall\": el operador xp de Berry–Keating")
    linea("  (cierre conceptual; NO es un intento de RH)")
    linea("=" * 72)

    # ------------------------------------------------------------------
    # A.1 — Silueta: el conteo semiclásico de xp = densidad suave de RvM
    # ------------------------------------------------------------------
    linea("\n## A.1 — xp captura la SILUETA (densidad suave)")
    Es = np.array([50.0, 200.0, 1000.0, 2500.0])
    linea("| E | xp semiclásico N(E) | Riemann–von Mangoldt suave | leading (E/2π)(logE−1) |")
    linea("|---|---------------------|----------------------------|------------------------|")
    for E in Es:
        linea(f"| {E:.0f} | {f(xp.conteo_semiclasico_xp(E))} | {f(spacing.N_suave(E))} | "
              f"{f(float(xp.conteo_leading_order(E)))} |")
    dif = float(np.max(np.abs(xp.conteo_semiclasico_xp(Es) - spacing.N_suave(Es))))
    linea(f"\n- diferencia xp-semiclásico vs RvM-suave: {dif:.2e} (son la MISMA fórmula).")
    linea("- Esta es exactamente la densidad que usa el unfolding del repo "
          "(`spacing.N_suave`). De ahí que xp sea el sospechoso natural.")

    # ------------------------------------------------------------------
    # A.2 — Rostro: el generador de dilatación discretizado da un peine
    # ------------------------------------------------------------------
    linea("\n## A.2 — el generador de dilatación discretizado: un PEINE uniforme")
    E_xp = xp.espectro_xp(N_GRID, L_CAJA)
    niveles = xp.niveles_positivos(E_xp, recorte=0.4)
    paso_teorico = 2 * np.pi / L_CAJA
    s_xp = xp.espaciados_peine(niveles)
    linea(f"- caja u∈[0,{L_CAJA}], grid {N_GRID} puntos; H = −i d/du (momento en u).")
    linea(f"- primeros niveles positivos: {np.round(niveles[:6], 4).tolist()}")
    linea(f"- paso del peine ≈ {np.diff(niveles).mean():.4f}  (teórico 2π/L = {paso_teorico:.4f})")
    linea(f"- desviación estándar de los espaciados unfolded: {s_xp.std():.5f} "
          "(≈0 ⇒ peine rígido, NO los t_n).")

    # ------------------------------------------------------------------
    # A.3 — comparación con los ceros y con GUE
    # ------------------------------------------------------------------
    gammas = zeros.cargar_ceros(N_ZEROS)
    s_zeta = spacing.unfold_densidad_local(gammas)
    s_gue = gue.espaciados_gue(2000, rng, n_matrices=10, recorte=0.1)
    m_xp, m_gue, m_zeta = (metrics.momentos(s) for s in (s_xp, s_gue, s_zeta))
    linea("\n## A.3 — niveles y estadística: xp vs ζ vs GUE")
    linea("| espectro | var espaciados | P(s<0.5) | repulsión (P(s<0.25)) |")
    linea("|----------|----------------|----------|------------------------|")
    linea(f"| xp (peine) | {f(m_xp['var'])} | {f(m_xp['frac_menor_0.5'])} | "
          f"{f(float((s_xp < 0.25).mean()))} |")
    linea(f"| GUE | {f(m_gue['var'])} | {f(m_gue['frac_menor_0.5'])} | "
          f"{f(float((s_gue < 0.25).mean()))} |")
    linea(f"| ζ (ceros) | {f(m_zeta['var'])} | {f(m_zeta['frac_menor_0.5'])} | "
          f"{f(float((s_zeta < 0.25).mean()))} |")
    linea("\n- El peine xp tiene varianza de espaciados ≈ 0 (rígido), mientras GUE y los "
          "ceros comparten la distribución de Wigner con repulsión. xp es dinámicamente "
          "demasiado pobre para generar estadística GUE.")

    # ------------------------------------------------------------------
    # Figuras
    # ------------------------------------------------------------------
    figs = [
        viz.silueta_conteo(gammas, f"{FIGDIR}/A1_silueta_conteo.png"),
        viz.peine_xp(niveles, L_CAJA, f"{FIGDIR}/A2_peine.png"),
        viz.niveles_vs_zeros(gammas, f"{FIGDIR}/A3_niveles_vs_zeros.png"),
        viz.espaciados_xp(s_xp, s_gue, s_zeta, f"{FIGDIR}/A3_espaciados.png"),
    ]
    linea("\n## Figuras")
    for fg in figs:
        linea(f"- {fg}")

    # ------------------------------------------------------------------
    # Veredicto A + tesis
    # ------------------------------------------------------------------
    linea("\n## Veredicto A")
    linea("xp reproduce la densidad media (silueta) EXACTAMENTE, pero su espectro discreto "
          "es un peine uniforme: no da los niveles individuales t_n, ni la repulsión tipo "
          "GUE, ni las correcciones aritméticas. **Toca la silueta, no el rostro.**")

    linea("\n## Tesis del reporte")
    linea("El operador buscado tendría que unificar TRES capas:")
    linea("  1. la densidad suave  → la da xp (silueta);")
    linea("  2. la estadística GUE → la dan las matrices aleatorias (Montgomery–Odlyzko), "
          "pero como modelo estadístico, no como un operador con esos autovalores;")
    linea("  3. la huella aritmética de los primos → vive en la fórmula explícita, no en xp.")
    linea("Los modelos existentes capturan PIEZAS, no el objeto. El programa de Hilbert–Pólya "
          "no está bloqueado por falta de cómputo, sino porque el objeto/geometría que "
          "tendría esas tres capas a la vez no está definido. El código ilumina el muro; "
          "no lo atraviesa. (Lenguaje prudente: esto es expositivo, no una contribución a RH.)")
    linea("\n" + "=" * 72)

    out = Path(__file__).resolve().parents[1] / "docs" / "riemann_missing_operator.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    cab = ["# El operador que falta — por qué xp \"casi\" funciona (Touching the Wall)", "",
           "_Demostración expositiva del Riemann Statistical Lab. NO es un intento de RH._",
           f"_Generado por `scripts/run_xp_demo.py` (grid {N_GRID}, L={L_CAJA}, semilla {SEED})._",
           ""]
    out.write_text("\n".join(cab + R), encoding="utf-8")
    print(f"\nReporte escrito en: {out}")


if __name__ == "__main__":
    main()
