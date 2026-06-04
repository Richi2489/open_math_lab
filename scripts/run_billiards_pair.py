"""
run_billiards_pair.py — el par controlado: integrable vs caótico (BGS)
======================================================================

Genera con código la conexión caos → repulsión de niveles (Bohigas–Giannoni–Schmit) que
el programa espectral de Riemann tomó prestada. MISMO solver, MISMA maquinaria estadística
(reusada del lab de Riemann); sólo cambia la geometría de la mesa:

  - rectángulo de lados inconmensurables (integrable)  → Poisson;
  - cuarto de estadio de Bunimovich (caótico, desimetrizado) → GOE.

Escribe docs/billiards_chaos_vs_integrable.md y figuras en docs/figures/riemann_billiards/.
Cachea el espectro del estadio en data/billiards/ (el FEM/FD es lo único caro).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, ValueError):  # pragma: no cover
    pass

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ))

from billiards import solver, viz, weyl  # noqa: E402
from riemann import gue, metrics  # noqa: E402

FIGDIR = "docs/figures/riemann_billiards"
# rectángulo: razón áurea. NO basta "inconmensurable" (a/b irracional): se necesita el
# CUADRADO de la razón irracional, o aparecen degeneraciones aritméticas (b=√2 da b²/a²=2
# racional → no-Poisson). φ²=φ+1 es irracional y "mal aproximable" → Poisson limpio.
A_REC, B_REC = 1.0, float((1.0 + np.sqrt(5.0)) / 2.0)
A_EST, R_EST = 1.0, 1.0                    # cuarto de estadio


def f(x, n=4):
    return f"{x:.{n}f}"


def estadio_eigs(a, r, h, n_modos):
    """Espectro del cuarto de estadio, cacheado (el solver es lo único caro)."""
    d = RAIZ / "data" / "billiards"
    npy, meta = d / "cuarto_estadio_eigs.npy", d / "cuarto_estadio_meta.json"
    params = {"a": a, "r": r, "h": h, "n_modos": n_modos}
    if npy.exists() and meta.exists() and json.loads(meta.read_text()) == params:
        return np.load(npy)
    print(f"[ resolviendo estadio (h={h}, {n_modos} modos)… cacheando ]")
    E = solver.resolver_dirichlet(
        lambda X, Y: solver.cuarto_estadio_dentro(X, Y, a, r), 0, a + r, 0, r, h, n_modos)
    d.mkdir(parents=True, exist_ok=True)
    np.save(npy, E)
    meta.write_text(json.dumps(params), encoding="utf-8")
    return E


def espaciados_norm(E, area, perim, descartar):
    """Espaciados unfolded por Weyl, descartando los modos bajos, normalizados a media 1.
    Devuelve (espaciados, media_pre_normalizacion)."""
    E = np.sort(np.asarray(E, dtype=np.float64))[descartar:]
    s = weyl.espaciados(E, area, perim)
    media = float(s.mean())
    return s / media, media


def stats(s, s_goe, s_pois):
    m = metrics.momentos(s)
    return {
        "n": len(s), "var": m["var"], "p05": m["frac_menor_0.5"],
        "p01": float((s < 0.1).mean()),
        "ks_goe": metrics.distancias(s, s_goe)["ks"],
        "wass_goe": metrics.distancias(s, s_goe)["wasserstein"],
        "ks_pois": metrics.distancias(s, s_pois)["ks"],
        "wass_pois": metrics.distancias(s, s_pois)["wasserstein"],
    }


def main():
    ap = argparse.ArgumentParser(description="Par controlado de billares (BGS).")
    ap.add_argument("--h", type=float, default=0.015, help="paso del grid del estadio.")
    ap.add_argument("--modos", type=int, default=700, help="modos del estadio.")
    ap.add_argument("--descartar", type=int, default=100, help="modos bajos descartados.")
    ap.add_argument("--n-rect", type=int, default=2000, help="modos del rectángulo (analítico).")
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--rapido", action="store_true")
    args = ap.parse_args()
    if args.rapido:
        args.h, args.modos, args.descartar, args.n_rect = 0.03, 200, 30, 800

    rng = np.random.default_rng(args.seed)
    R = []

    def linea(s=""):
        print(s)
        R.append(s)

    linea("=" * 72)
    linea("  Billares cuánticos — par controlado: integrable vs caótico (BGS)")
    linea("=" * 72)

    # -- Reglas pre-registradas --------------------------------------------
    linea("\n## Reglas pre-registradas (antes de calcular)")
    linea("- **Baseline del estadio = GOE** (hay simetría de inversión temporal), NO GUE: "
          "repulsión lineal P(s)∼s, var≈0.27. Comparar contra GUE mostraría una desviación "
          "ESPERADA, no una anomalía.")
    linea("- **Desimetrizar:** se trabaja en un solo sector (cuarto de estadio, Dirichlet en "
          "todo el borde, ejes de simetría incluidos). Mezclar sectores destruye la repulsión "
          "y finge Poisson — es el confounder principal del experimento.")
    linea("- **Unfolding por Weyl con término de perímetro:** N(E)≈(A/4π)E−(P/4π)√E, no sólo "
          "el área.")
    linea("- **Tamaño finito:** los modos bajos se desvían; el match mejora alto en el "
          f"espectro (se descartan los primeros {args.descartar}).")

    # -- Dominios y espectros ----------------------------------------------
    area_r, perim_r = weyl.geom_rectangulo(A_REC, B_REC)
    area_e, perim_e = weyl.geom_cuarto_estadio(A_EST, R_EST)
    E_rect = solver.eigs_rectangulo(A_REC, B_REC, args.n_rect)
    E_estad = estadio_eigs(A_EST, R_EST, args.h, args.modos)

    s_rect, mr = espaciados_norm(E_rect, area_r, perim_r, args.descartar)
    s_estad, me = espaciados_norm(E_estad, area_e, perim_e, args.descartar)

    linea("\n## Dominios (mismo solver de Helmholtz–Dirichlet)")
    linea(f"- rectángulo {A_REC}×{B_REC:.4f} (razón áurea, cuadrado irracional → evita "
          f"degeneraciones aritméticas): {len(E_rect)} modos analíticos; A={area_r:.3f}, "
          f"P={perim_r:.3f}.")
    linea(f"- cuarto de estadio a={A_EST}, r={R_EST} (FD h={args.h}): {len(E_estad)} modos; "
          f"A={area_e:.3f}, P={perim_e:.3f}; E_max={E_estad[-1]:.0f}.")
    linea(f"- media de espaciado pre-normalización (Weyl): rectángulo {mr:.3f}, estadio "
          f"{me:.3f}. La leve desviación del estadio refleja la frontera escalonada del FD; "
          "se normaliza a media 1 (paso estándar de unfolding).")

    # -- Baselines ---------------------------------------------------------
    s_goe = gue.espaciados_goe(300, rng, n_matrices=30, recorte=0.1)
    s_pois = gue.espaciados_poisson(len(s_goe), rng)

    # -- Estadística -------------------------------------------------------
    st_r = stats(s_rect, s_goe, s_pois)
    st_e = stats(s_estad, s_goe, s_pois)
    linea("\n## Estadística de espaciados (misma maquinaria del lab de Riemann)")
    linea("| dominio | n | var | P(s<0.5) | P(s<0.1) | KS vs Poisson | KS vs GOE |")
    linea("|---------|---|-----|----------|----------|---------------|-----------|")
    linea(f"| rectángulo (integrable) | {st_r['n']} | {f(st_r['var'])} | {f(st_r['p05'])} | "
          f"{f(st_r['p01'])} | {f(st_r['ks_pois'])} | {f(st_r['ks_goe'])} |")
    linea(f"| estadio (caótico) | {st_e['n']} | {f(st_e['var'])} | {f(st_e['p05'])} | "
          f"{f(st_e['p01'])} | {f(st_e['ks_pois'])} | {f(st_e['ks_goe'])} |")
    linea(f"| — Poisson (ref.) | {len(s_pois)} | {f(s_pois.var())} | "
          f"{f(float((s_pois < 0.5).mean()))} | {f(float((s_pois < 0.1).mean()))} | — | — |")
    linea(f"| — GOE (ref.) | {len(s_goe)} | {f(s_goe.var())} | "
          f"{f(float((s_goe < 0.5).mean()))} | {f(float((s_goe < 0.1).mean()))} | — | — |")

    linea("\nLectura: el rectángulo está más cerca de Poisson que de GOE; el estadio, al "
          "revés. P(s<0.1) ≈ 0 en el estadio (hoyo de repulsión); grande en el rectángulo.")

    # -- Figuras -----------------------------------------------------------
    figs = [
        viz.par_espaciados(s_rect, s_estad, f"{FIGDIR}/par_espaciados.png"),
        viz.cdf_par(s_rect, s_estad, s_goe, s_pois, f"{FIGDIR}/cdf_par.png"),
        viz.conteo_weyl(E_estad, area_e, perim_e, f"{FIGDIR}/weyl_estadio.png", "cuarto de estadio"),
    ]
    linea("\n## Figuras")
    for fg in figs:
        linea(f"- {fg}")

    # -- Veredicto ---------------------------------------------------------
    rect_poisson = st_r["ks_pois"] < st_r["ks_goe"] and st_r["var"] > 0.7
    estadio_goe = st_e["ks_goe"] < st_e["ks_pois"] and st_e["var"] < 0.45
    repulsion_por_geometria = st_e["p05"] < st_r["p05"] and st_e["var"] < st_r["var"]
    linea("\n## Veredicto (lenguaje prudente)")
    linea(f"- rectángulo ≈ Poisson (no GOE): {'SÍ' if rect_poisson else 'NO'} "
          f"(KS vs Poisson {f(st_r['ks_pois'])} < KS vs GOE {f(st_r['ks_goe'])}; var≈1)")
    linea(f"- estadio ≈ GOE (no Poisson): {'SÍ' if estadio_goe else 'NO'} "
          f"(KS vs GOE {f(st_e['ks_goe'])} < KS vs Poisson {f(st_e['ks_pois'])}; var≈0.27)")
    linea(f"- la repulsión aparece SÓLO al volver caótica la mesa: "
          f"{'SÍ' if repulsion_por_geometria else 'NO'} "
          f"(P(s<0.5): {f(st_r['p05'])} → {f(st_e['p05'])}; var: {f(st_r['var'])} → {f(st_e['var'])})")
    if rect_poisson and estadio_goe and repulsion_por_geometria:
        linea("\n**SE REPRODUCE BGS.** Con el MISMO solver y la MISMA maquinaria estadística, "
              "cambiar la geometría de integrable a caótica hace aparecer la repulsión de "
              "niveles tipo GOE; el caso integrable se queda en Poisson. La conexión "
              "caos → repulsión que el lab de Riemann tomó prestada queda GENERADA aquí desde "
              "cero. Prudencia: es una demostración de Bohigas–Giannoni–Schmit, no teoría nueva, "
              "y el match es de tamaño finito (mejora alto en el espectro).")
    else:
        linea("\n**Revisar** alguna condición (desimetrización, unfolding o nº de modos) antes "
              "de leer el resultado.")
    linea("\n" + "=" * 72)

    out = RAIZ / "docs" / "billiards_chaos_vs_integrable.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    cab = ["# Billares cuánticos — caos vs integrabilidad (Bohigas–Giannoni–Schmit)", "",
           "_Demostración del open_math_lab. NO es teoría nueva._",
           f"_Generado por `scripts/run_billiards_pair.py` (estadio h={args.h}, "
           f"{args.modos} modos, semilla {args.seed})._", ""]
    out.write_text("\n".join(cab + R), encoding="utf-8")
    print(f"\nReporte escrito en: {out}")


if __name__ == "__main__":
    main()
