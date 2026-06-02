"""
run_odlyzko_highT.py — confirmación a alta altura con ceros de Odlyzko (iteración 2)
===================================================================================

Confirma (o no) el cierre de la réplica Montgomery–Odlyzko subiendo la altura de ~10³
(mpmath, 2000 ceros) a ~10⁶ (Odlyzko, ~2M ceros). NADA de proof-hunting.

TRAMPA PRE-REGISTRADA: con ~2M de espaciados el KS se vuelve hipersensible y casi seguro
marcará p<0.05. ESO NO ES ANOMALÍA: es la corrección conocida de altura finita al GUE. La
regla manda leer la TENDENCIA (¿la desviación encoge respecto de baja altura?) y comparar
contra GUE de tamaño finito, NO cazar el p-valor del KS.

Uso:
    python scripts/run_odlyzko_highT.py --path data/odlyzko/zeros6 --n 2000000 \
        --windows 5 --gue-dim 2000 --gue-samples 40 --seed 2026
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

from riemann import gue, metrics, odlyzko, spacing, viz, zeros  # noqa: E402

FIGDIR = "docs/figures/riemann_odlyzko"


def f(x, n=4):
    return f"{x:.{n}f}"


def unfold_zeta(gammas, metodo):
    if metodo == "theta_exact":
        return np.diff(spacing.theta_unfold(gammas))
    return spacing.unfold_densidad_local(gammas)


def bundle(s, s_gue, s_pois):
    m = metrics.momentos(s)
    dg = metrics.distancias(s, s_gue)
    dp = metrics.distancias(s, s_pois)
    return {
        "n": len(s), "media": m["media"], "var": m["var"],
        "p05": m["frac_menor_0.5"], "p025": float((np.asarray(s) < 0.25).mean()),
        "ks_gue": dg["ks"], "ks_gue_p": dg["ks_p"], "wass_gue": dg["wasserstein"],
        "ks_pois": dp["ks"], "wass_pois": dp["wasserstein"],
    }


def main():
    ap = argparse.ArgumentParser(description="Confirmación alta-T con Odlyzko.")
    ap.add_argument("--path", default="data/odlyzko/zeros6")
    ap.add_argument("--n", type=int, default=2_000_000)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--windows", type=int, default=5)
    ap.add_argument("--gue-dim", type=int, default=2000)
    ap.add_argument("--gue-samples", type=int, default=40)
    ap.add_argument("--unfold", choices=["local_density", "theta_exact"], default="local_density")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--out", default="docs/odlyzko_highT_report.md")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    R = []

    def linea(s=""):
        print(s)
        R.append(s)

    linea("=" * 72)
    linea("  Riemann Lab — confirmación alta-T con ceros de Odlyzko (iteración 2)")
    linea("=" * 72)

    # -- A. Regla pre-registrada -------------------------------------------
    linea("\n## A. Regla pre-registrada (con la trampa del KS)")
    linea("Comparar SIEMPRE contra GUE de tamaño finito simulado, no contra la fórmula "
          "asintótica. **Trampa del KS:** con ~2M de muestras el KS es hipersensible y casi "
          "seguro dará p<0.05; eso NO es anomalía, es la corrección conocida de altura finita "
          "al GUE. El veredicto se lee por la TENDENCIA (¿la desviación encoge vs baja altura "
          "y por ventanas?), por la separación clara frente a Poisson, y por si el residuo "
          "cabe en la variabilidad del GUE finito — NUNCA por el p-valor del KS.")

    # -- B. Dataset y validación -------------------------------------------
    print("\n[ cargando Odlyzko... ]")
    g_all = odlyzko.cargar_odlyzko(args.path)
    odlyzko.validar_zeros(g_all)
    g = odlyzko.tomar_ventana(g_all, args.start, n=args.n)
    odlyzko.validar_zeros(g)
    res = odlyzko.resumir_dataset(g)
    linea("\n## B. Dataset y validaciones")
    linea(f"- archivo: `{args.path}` | total en archivo: {len(g_all):,} ceros")
    linea(f"- muestra usada: N = {res['n']:,} (start={args.start})")
    linea(f"- altura T: {res['T_min']:.3f} … {res['T_max']:.1f}")
    linea(f"- validación: 1D numérico, estrictamente ascendente, sin duplicados/NaN/inf ✓")

    # -- C. Unfolding -------------------------------------------------------
    linea("\n## C. Método de unfolding")
    linea(f"- usado: **{args.unfold}**"
          + ("  (s_n = (γ_{n+1}−γ_n)·log(γ_n/2π)/2π)" if args.unfold == "local_density"
             else "  (θ(γ)/π+1 exacto; solo control de muestra chica)"))
    s_zeta = unfold_zeta(g, args.unfold)

    # Baselines GUE y Poisson (tamaño finito).
    print("[ simulando GUE finito... ]")
    s_gue = gue.espaciados_gue(args.gue_dim, rng, n_matrices=args.gue_samples, recorte=0.1)
    s_pois = gue.espaciados_poisson(len(s_gue), rng)
    var_gue, p05_gue = float(s_gue.var()), float((s_gue < 0.5).mean())

    # -- D. Sanity checks ---------------------------------------------------
    linea("\n## D. Sanity checks")
    linea(f"- media del espaciado unfolded (debe ≈ 1): {s_zeta.mean():.5f}")
    g_mp = zeros.cargar_ceros(2000)
    g_od_2000 = odlyzko.tomar_ventana(g_all, 0, n=2000)
    dif = float(np.max(np.abs(g_mp - g_od_2000)))
    linea(f"- cross-check Odlyzko vs mpmath (primeros 2000 γ): máx |Δ| = {dif:.2e} "
          f"(precisión Odlyzko ~4e-9) {'✓' if dif < 1e-6 else '⚠'}")
    linea(f"- GUE finito ({args.gue_dim}×{args.gue_dim}, {len(s_gue):,} espac.): "
          f"var={var_gue:.4f}, P(s<0.5)={p05_gue:.4f}")

    # -- E. ζ alta-T vs GUE / Poisson / Wigner ------------------------------
    print("[ métricas alta-T... ]")
    b = bundle(s_zeta, s_gue, s_pois)
    linea("\n## E. ζ alta-T vs GUE, Poisson, Wigner")
    linea(f"- ζ (N={b['n']:,}): media={f(b['media'])}, var={f(b['var'])}, "
          f"P(s<0.5)={f(b['p05'])}, P(s<0.25)={f(b['p025'])}")
    linea("| comparación | KS | KS p-valor | Wasserstein |")
    linea("|-------------|----|-----------|-------------|")
    linea(f"| ζ vs GUE | {f(b['ks_gue'])} | {b['ks_gue_p']:.2e} | {f(b['wass_gue'])} |")
    linea(f"| ζ vs Poisson | {f(b['ks_pois'])} | {metrics.distancias(s_zeta, s_pois)['ks_p']:.2e} | {f(b['wass_pois'])} |")
    linea(f"- KS p-valor ζ-GUE = {b['ks_gue_p']:.2e}: con N={b['n']:,} esto es la TRAMPA "
          "esperada (hipersensibilidad), no una anomalía. Lo informativo es el TAMAÑO de la "
          f"desviación (KS={f(b['ks_gue'])}, Wasserstein={f(b['wass_gue'])}) y su tendencia.")

    # -- F. Baja altura vs Odlyzko ------------------------------------------
    s_mp = unfold_zeta(g_mp, args.unfold)
    s_od2000 = unfold_zeta(g_od_2000, args.unfold)
    linea("\n## F. Baja altura vs Odlyzko alta-T")
    linea("| dataset | N spac | T_min | T_max | var | P(s<0.5) | KS vs GUE | Wass vs GUE | KS vs Poisson |")
    linea("|---------|--------|-------|-------|-----|----------|-----------|-------------|---------------|")
    filas_F = [
        ("mpmath N=2000", g_mp, s_mp),
        ("Odlyzko primeros 2000", g_od_2000, s_od2000),
        (f"Odlyzko N={args.n:,}", g, s_zeta),
    ]
    bundles_F = {}
    for nombre, gg, ss in filas_F:
        bi = bundle(ss, s_gue, s_pois)
        bundles_F[nombre] = bi
        linea(f"| {nombre} | {bi['n']:,} | {gg[0]:.0f} | {gg[-1]:.0f} | {f(bi['var'])} | "
              f"{f(bi['p05'])} | {f(bi['ks_gue'])} | {f(bi['wass_gue'])} | {f(bi['ks_pois'])} |")
    linea(f"\n(GUE finito de referencia: var={var_gue:.4f}, P(s<0.5)={p05_gue:.4f}. "
          "Leer cómo |var−var_GUE| y |P(s<0.5)−P_GUE| cambian al subir la altura.)")

    # -- G. Estabilidad por ventanas ----------------------------------------
    print("[ ventanas... ]")
    ventanas = odlyzko.dividir_en_ventanas(g, args.windows)
    linea("\n## G. Estabilidad por ventanas no solapadas (altura creciente)")
    linea("| ventana | T_min | T_max | N spac | var | P(s<0.5) | |Δvar| | KS vs GUE | Wass vs GUE |")
    linea("|---------|-------|-------|--------|-----|----------|--------|-----------|-------------|")
    cT, ks_w, wass_w, dvar_w, dp05_w = [], [], [], [], []
    for i, vg in enumerate(ventanas, 1):
        sv = unfold_zeta(vg, args.unfold)
        bi = bundle(sv, s_gue, s_pois)
        dvar = abs(bi["var"] - var_gue)
        cT.append(0.5 * (vg[0] + vg[-1])); ks_w.append(bi["ks_gue"])
        wass_w.append(bi["wass_gue"]); dvar_w.append(dvar)
        dp05_w.append(abs(bi["p05"] - p05_gue))
        linea(f"| W{i} | {vg[0]:.0f} | {vg[-1]:.0f} | {bi['n']:,} | {f(bi['var'])} | "
              f"{f(bi['p05'])} | {f(dvar)} | {f(bi['ks_gue'])} | {f(bi['wass_gue'])} |")
    pend_dvar = float(np.polyfit(cT, dvar_w, 1)[0])
    pend_wass = float(np.polyfit(cT, wass_w, 1)[0])
    linea(f"\n- tendencia |Δvar| con altura: pendiente {pend_dvar:+.2e} "
          f"(W1={dvar_w[0]:.4f} → W{args.windows}={dvar_w[-1]:.4f})")
    linea(f"- tendencia Wasserstein con altura: pendiente {pend_wass:+.2e} "
          f"(W1={wass_w[0]:.4f} → W{args.windows}={wass_w[-1]:.4f})")

    # -- Figuras ------------------------------------------------------------
    print("[ figuras... ]")
    figs = [
        viz.espaciados_nn(s_zeta, s_gue, s_pois, f"{FIGDIR}/hist_spacing.png"),
        viz.cdf_espaciados(s_zeta, s_gue, s_pois, f"{FIGDIR}/cdf.png"),
        viz.cdf_diferencia(s_zeta, s_gue, f"{FIGDIR}/cdf_diferencia.png"),
        viz.metricas_por_ventana(cT, ks_w, wass_w, f"{FIGDIR}/por_ventana.png"),
        viz.comparar_datasets(
            ["mpmath\n2000", "Odlyzko\n2000", f"Odlyzko\n{args.n//1000}k"],
            [bundles_F[n]["wass_gue"] for n, _, _ in filas_F],
            f"{FIGDIR}/baja_vs_alta.png"),
    ]
    centros, r2 = metrics.correlacion_pares(np.cumsum(s_zeta), r_max=3.0, dr=0.1)
    figs.append(viz.pares_r2(centros, r2, f"{FIGDIR}/r2.png"))
    linea("\n## Figuras")
    for fg in figs:
        linea(f"- {fg}")

    # -- H. Limitaciones ----------------------------------------------------
    linea("\n## H. Limitaciones")
    linea("- Altura ~10⁶, no ~10²⁰: los datasets de muy alta altura (offsets) usan otro "
          "formato/parser; quedan para una iteración futura.")
    linea("- GUE finito de dim fija como referencia (cuasi-asintótica para el spacing NN); "
          "no se modela la correspondencia exacta dim↔altura.")
    linea("- KS hipersensible a N grande (trampa pre-registrada); por eso el veredicto se "
          "apoya en tamaño de desviación y tendencia, no en el p-valor.")
    linea("- R₂ con unfolding local (cumsum de espaciados), suficiente para la firma de "
          "repulsión, no para correlaciones de largo alcance finas.")

    # -- I. Veredicto -------------------------------------------------------
    lo = bundles_F["mpmath N=2000"]
    hi = bundles_F[f"Odlyzko N={args.n:,}"]
    dvar_lo, dvar_hi = abs(lo["var"] - var_gue), abs(hi["var"] - var_gue)
    poisson_descartado = hi["ks_pois"] > 5 * hi["ks_gue"] and hi["wass_pois"] > 3 * hi["wass_gue"]
    gue_baseline_ok = abs(hi["var"] - var_gue) < 0.03 and abs(hi["p05"] - p05_gue) < 0.03
    desviacion_encoge = dvar_hi < dvar_lo
    ventanas_estables = pend_dvar <= 1e-8  # no crece con la altura
    linea("\n## I. Veredicto (lenguaje prudente)")
    linea(f"- ¿la desviación a GUE ENCOGE vs baja altura? {'SÍ' if desviacion_encoge else 'NO'} "
          f"(|Δvar|: {dvar_lo:.4f} a alturas ~10³ → {dvar_hi:.4f} a ~10⁶)")
    linea(f"- ¿Poisson claramente descartado? {'SÍ' if poisson_descartado else 'NO'} "
          f"(KS ζ-Poisson={f(hi['ks_pois'])} ≫ KS ζ-GUE={f(hi['ks_gue'])})")
    linea(f"- ¿GUE sigue siendo el baseline correcto? {'SÍ' if gue_baseline_ok else 'NO'} "
          f"(var y P(s<0.5) de ζ ≈ GUE finito)")
    linea(f"- ¿desviación estable/decreciente por ventanas (no residuo creciente)? "
          f"{'SÍ' if ventanas_estables else 'NO'}")
    linea(f"- KS p-valor ζ-GUE = {hi['ks_gue_p']:.1e}: trampa esperada, NO se interpreta como anomalía.")

    cierre = poisson_descartado and gue_baseline_ok and desviacion_encoge and ventanas_estables
    if cierre:
        linea("\n**SE CONFIRMA EL CIERRE DE LA RÉPLICA.** A altura ~10⁶ los ceros de ζ están "
              "mucho más cerca del GUE finito que de Poisson, y la pequeña desviación de baja "
              "altura ENCOGE al subir la altura (consistente con la corrección conocida de "
              "tamaño/altura finita). No hay evidencia de anomalía estructural; el p-valor del "
              "KS es la trampa de muestra grande, no señal. La repulsión tipo GUE de los gaps "
              "unfolded de ζ queda reproducida limpiamente. Prudencia: es una réplica empírica "
              "de un resultado conocido (Montgomery–Odlyzko), no evidencia sobre RH.")
    else:
        linea("\n**NO se confirma limpiamente** alguna condición. Revisar la sección "
              "correspondiente antes de interpretar; la regla prohíbe leer una anomalía sin "
              "que sobreviva a tendencia, ventanas y baseline.")
    linea("\n" + "=" * 72)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cab = ["# Confirmación alta-T con ceros de Odlyzko — iteración 2", "",
           f"_Generado por `scripts/run_odlyzko_highT.py` (N={args.n}, unfold={args.unfold}, "
           f"GUE {args.gue_dim}×{args.gue_dim}×{args.gue_samples}, semilla {args.seed})._", ""]
    out.write_text("\n".join(cab + R), encoding="utf-8")
    print(f"\nReporte escrito en: {out}")


if __name__ == "__main__":
    main()
