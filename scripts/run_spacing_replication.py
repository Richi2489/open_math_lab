"""
run_spacing_replication.py — réplica Montgomery–Odlyzko (iteración 1)
====================================================================

¿Los gaps unfolded de los ceros de ζ exhiben repulsión tipo GUE? Sanity check
reproducible, contrastado contra GUE finito simulado y contra Poisson. NADA de claims
grandiosos; lenguaje prudente.

Imprime la REGLA anti-autoengaño ANTES de calcular y escribe
docs/spacing_replication_report.md con tablas y figuras.

Uso:
    python scripts/run_spacing_replication.py                 # N=2000 ceros
    python scripts/run_spacing_replication.py --n-ceros 5000
    python scripts/run_spacing_replication.py --rapido        # N chico para fumar el pipeline
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, ValueError):  # pragma: no cover
    pass

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from riemann import gue, metrics, spacing, viz, zeros  # noqa: E402


def fmt(x, n=4):
    return f"{x:.{n}f}"


def main():
    ap = argparse.ArgumentParser(description="Réplica Montgomery–Odlyzko (iteración 1).")
    ap.add_argument("--n-ceros", type=int, default=2000)
    ap.add_argument("--n-gue", type=int, default=0, help="dim. de la matriz GUE (0 = igual a N).")
    ap.add_argument("--ventanas", type=int, default=5)
    ap.add_argument("--semilla", type=int, default=20260602)
    ap.add_argument("--rapido", action="store_true")
    args = ap.parse_args()

    n_ceros = 400 if args.rapido else args.n_ceros
    n_gue = args.n_gue or min(n_ceros, 2000)
    rng = np.random.default_rng(args.semilla)
    R = []

    def linea(s=""):
        print(s)
        R.append(s)

    linea("=" * 70)
    linea("  Riemann Statistical Lab — réplica Montgomery–Odlyzko (iteración 1)")
    linea("=" * 70)

    # -- Regla anti-autoengaño (ANTES de calcular) --------------------------
    linea("\n## Regla anti-autoengaño (PRE-REGISTRADA, antes de ver resultados)")
    linea("No se declarará una desviación de GUE si desaparece al cambiar (a) la ventana de "
          "altura T, (b) el método de unfolding, (c) el tamaño de muestra, o (d) el baseline. "
          "Toda diferencia contra GUE se compara contra GUE de tamaño finito simulado (mismo N), "
          "NO contra la fórmula asintótica ideal. El 'éxito' de la iteración 1 es REPRODUCIR "
          "limpiamente la repulsión tipo GUE, no encontrar una anomalía.")

    # -- Parámetros ---------------------------------------------------------
    linea("\n## Parámetros")
    linea(f"- ceros de ζ: N = {n_ceros:,} (mpmath, cacheados en data/)")
    linea(f"- GUE finito: matriz {n_gue}×{n_gue}; baseline Poisson (gaps exp i.i.d.)")
    linea(f"- ventanas de altura: {args.ventanas} | semilla: {args.semilla}")

    # -- Ceros + unfolding --------------------------------------------------
    print("\n[ cargando ceros de ζ (puede tardar la primera vez)... ]")
    gammas = zeros.cargar_ceros(n_ceros, verbose=True)
    zeros.validar_ceros(gammas)
    s_zeta = spacing.gaps_vecino(gammas)
    w = spacing.unfold(gammas)
    linea(f"\n## Ceros y unfolding")
    linea(f"- γ_1 = {gammas[0]:.4f}, γ_2 = {gammas[1]:.4f}, γ_N = {gammas[-1]:.4f}")
    linea(f"- espaciados unfolded: {len(s_zeta):,} | media = {s_zeta.mean():.4f} "
          f"(debe ≈ 1 si el unfolding es correcto)")

    # -- Baselines ----------------------------------------------------------
    n_mat = max(1, int(np.ceil(3 * len(s_zeta) / (n_gue * 0.8))))
    s_gue = gue.espaciados_gue(n_gue, rng, n_matrices=n_mat, recorte=0.1)
    s_pois = gue.espaciados_poisson(len(s_gue), rng)

    # -- A. Momentos --------------------------------------------------------
    mz, mg, mp = (metrics.momentos(s) for s in (s_zeta, s_gue, s_pois))
    linea("\n## A. Momentos de los espaciados")
    linea("| dist. | media | var | asimetría | P(s<0.5) |")
    linea("|-------|-------|-----|-----------|----------|")
    for nombre, m in [("ζ", mz), ("GUE", mg), ("Poisson", mp)]:
        linea(f"| {nombre} | {fmt(m['media'])} | {fmt(m['var'])} | {fmt(m['asimetria'])} | "
              f"{fmt(m['frac_menor_0.5'])} |")
    linea("P(s<0.5) chico = repulsión (pocos gaps diminutos). Poisson ~0.39; GUE mucho menor.")

    # -- B. Distancias vs GUE y vs Poisson ----------------------------------
    d_zg = metrics.distancias(s_zeta, s_gue)
    d_zp = metrics.distancias(s_zeta, s_pois)
    linea("\n## B. Distancias de la distribución de espaciados")
    linea("| comparación | KS | KS p-valor | Wasserstein |")
    linea("|-------------|----|-----------|-------------|")
    linea(f"| ζ vs GUE | {fmt(d_zg['ks'])} | {d_zg['ks_p']:.3e} | {fmt(d_zg['wasserstein'])} |")
    linea(f"| ζ vs Poisson | {fmt(d_zp['ks'])} | {d_zp['ks_p']:.3e} | {fmt(d_zp['wasserstein'])} |")
    linea("Repulsión tipo GUE ⇔ ζ MUCHO más cerca de GUE que de Poisson (KS y Wasserstein).")

    # -- C. Correlación de pares R_2 ----------------------------------------
    centros, r2 = metrics.correlacion_pares(w, r_max=3.0, dr=0.1)
    linea("\n## C. Correlación de pares R₂ (firma Montgomery–Odlyzko)")
    linea("| r | R₂ empírico (ζ) | GUE 1−(sinπr/πr)² |")
    linea("|---|-----------------|--------------------|")
    for r_obj in (0.2, 0.5, 1.0, 1.5, 2.0):
        i = int(np.argmin(np.abs(centros - r_obj)))
        linea(f"| {centros[i]:.1f} | {fmt(r2[i])} | {fmt(float(metrics.r2_teorico_gue(centros[i])))} |")
    linea("La firma es el HOYO de repulsión cerca de r→0 (R₂→0), no el valor plano 1 de Poisson.")

    # -- D. Confounder: unfolding INCORRECTO --------------------------------
    s_mal = spacing.gaps_normalizacion_global_INCORRECTA(gammas)
    m_mal = metrics.momentos(s_mal)
    d_mal = metrics.distancias(s_mal, s_gue)
    linea("\n## D. El confounder: normalización GLOBAL (INCORRECTA)")
    linea(f"- P(s<0.5) con unfolding correcto = {mz['frac_menor_0.5']:.4f} ; "
          f"con normalización global = {m_mal['frac_menor_0.5']:.4f}")
    linea(f"- var correcto = {mz['var']:.4f} ; var global-incorrecto = {m_mal['var']:.4f}  "
          "(la densidad cambia con la altura: dividir por la media global mezcla regímenes "
          "y distorsiona la distribución — el análogo del confounder de magnitud en Collatz).")

    # -- E. Estabilidad por altura ------------------------------------------
    estab = metrics.estabilidad_por_altura(gammas, args.ventanas, rng)
    linea("\n## E. Estabilidad por ventana de altura (efectos de tamaño finito)")
    linea("| ventana T | n espac. | KS(ζ vs GUE) | Wasserstein |")
    linea("|-----------|----------|--------------|-------------|")
    for e in estab:
        linea(f"| {e['T_min']:.0f}–{e['T_max']:.0f} | {e['n_spac']:,} | "
              f"{fmt(e['ks'])} | {fmt(e['wasserstein'])} |")
    ks_vals = [e["ks"] for e in estab]
    linea(f"\nKS por ventana: min {min(ks_vals):.3f}, max {max(ks_vals):.3f}. A baja altura "
          "el match suele ser peor (efecto de tamaño finito esperado, no anomalía).")

    # -- Figuras ------------------------------------------------------------
    print("\n[ guardando figuras... ]")
    figs = [
        viz.espaciados_nn(s_zeta, s_gue, s_pois),
        viz.cdf_espaciados(s_zeta, s_gue, s_pois),
        viz.pares_r2(centros, r2),
        viz.estabilidad_altura(estab),
    ]
    linea("\n## Figuras")
    for f in figs:
        linea(f"- {f}")

    # -- Veredicto ----------------------------------------------------------
    mas_cerca_gue = d_zg["ks"] < d_zp["ks"] and d_zg["wasserstein"] < d_zp["wasserstein"]
    repulsion = mz["frac_menor_0.5"] < 0.5 * mp["frac_menor_0.5"]
    r2_hoyo = r2[int(np.argmin(np.abs(centros - 0.2)))] < 0.6
    linea("\n## Veredicto (iteración 1) — lenguaje prudente")
    linea(f"- ζ más cerca de GUE que de Poisson (KS y Wasserstein): {'SÍ' if mas_cerca_gue else 'NO'}")
    linea(f"- repulsión a corta distancia (P(s<0.5) ≪ Poisson):     {'SÍ' if repulsion else 'NO'}")
    linea(f"- R₂ con hoyo de repulsión cerca de r→0:                {'SÍ' if r2_hoyo else 'NO'}")
    if mas_cerca_gue and repulsion and r2_hoyo:
        linea("\n**Repulsión tipo GUE reproducida** en las tres métricas (espaciados NN, "
              "distancias vs baselines, y correlación de pares). El residuo frente al GUE "
              "asintótico es compatible con tamaño finito / baja altura (ver sección E), no "
              "se interpreta como anomalía. Prudencia: es una réplica empírica de una señal "
              "conocida, no evidencia sobre RH.")
    else:
        linea("\n**Repulsión NO reproducida limpiamente** en alguna métrica. Antes de "
              "interpretar nada, revisar unfolding, tamaño de muestra y baseline (la regla "
              "pre-registrada exige que cualquier desviación sobreviva a esos cambios).")
    linea("\n" + "=" * 70)

    docs = Path(__file__).resolve().parents[1] / "docs"
    docs.mkdir(exist_ok=True)
    rep = docs / "spacing_replication_report.md"
    cab = ["# Réplica Montgomery–Odlyzko — iteración 1", "",
           f"_Generado por `scripts/run_spacing_replication.py` "
           f"(N={n_ceros}, semilla {args.semilla}{', rápido' if args.rapido else ''})._", ""]
    rep.write_text("\n".join(cab + R), encoding="utf-8")
    print(f"\nReporte escrito en: {rep}")


if __name__ == "__main__":
    main()
