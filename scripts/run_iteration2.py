"""
run_iteration2.py — ¿sobrevive el lag-1 positivo?
=================================================

Decide, con rigor falsable, si la autocorrelación de lag-1 de la secuencia de v's
del mapa acelerado es señal real de corto alcance o artefacto de muestra finita.

Hace, en orden:
  A. Sanity check (E[v], drift, distribución de v).
  B. ACF bruto agregado (lags 1..10).
  C. Corrección por sesgo finito -1/(L-1) y banda analítica.
  D. Prueba de PERMUTACIÓN (principal) por lag.
  E. Estratificación por longitud L (confound de supervivencia).
  F. Escalamiento por bits (¿lag-1 estable o se desvanece?).
  G. Veredicto explícito según la regla de decisión.

Escribe un reporte en docs/ y figuras en outputs/. Lenguaje estadístico prudente.

Uso:
    python scripts/run_iteration2.py             # corrida completa
    python scripts/run_iteration2.py --rapido    # versión chica para fumar el pipeline
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

import numpy as np

try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, ValueError):  # pragma: no cover
    pass

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz import acf, sampling, stats, viz  # noqa: E402

LAGS = list(range(1, 11))
BANDAS_BITS = [(20, 24), (30, 34), (40, 44), (50, 54), (60, 64)]
BANDA_PRINCIPAL = (40, 44)
MIN_TRAYECTORIAS_BUCKET = 200  # mínimo para correr permutación en un bucket


# ---------------------------------------------------------------------------
# Utilidades de formato
# ---------------------------------------------------------------------------
def fmt(x, n=4):
    return f"{x:+.{n}f}"


def construir_config(rapido: bool, semilla: int) -> dict:
    if rapido:
        return dict(
            semilla=semilla,
            B=300,
            n_eff_principal=40_000,
            n_eff_bits=30_000,
            n_drift=100_000,
        )
    return dict(
        semilla=semilla,
        B=1000,
        n_eff_principal=300_000,
        n_eff_bits=200_000,
        n_drift=1_000_000,
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="Iteración 2: ¿sobrevive el lag-1?")
    parser.add_argument("--rapido", action="store_true", help="versión chica y rápida.")
    parser.add_argument("--semilla", type=int, default=20260601)
    args = parser.parse_args()
    cfg = construir_config(args.rapido, args.semilla)

    rng_py = random.Random(cfg["semilla"])
    rng_np = np.random.default_rng(cfg["semilla"] + 1)
    R = []  # líneas del reporte markdown

    def linea(s=""):
        print(s)
        R.append(s)

    linea("=" * 70)
    linea("  open_math_lab — Collatz, iteración 2")
    linea("  ¿sobrevive el lag-1 positivo tras corrección Y permutación?")
    linea("=" * 70)

    # -- 1. Reproducibilidad ------------------------------------------------
    linea("\n## 1. Reproducibilidad / parámetros")
    linea(f"- semilla maestra: {cfg['semilla']}")
    linea(f"- permutaciones B: {cfg['B']}")
    linea(f"- N_eff objetivo (banda principal {BANDA_PRINCIPAL} bits): {cfg['n_eff_principal']:,}")
    linea(f"- N_eff objetivo por banda de bits: {cfg['n_eff_bits']:,}")
    linea(f"- lags analizados: {LAGS[0]}..{LAGS[-1]}")
    linea(f"- min_L (para soportar todos los lags): {LAGS[-1] + 1}")
    linea(f"- enteros: nativos de Python (sin riesgo de overflow)")

    # -- Muestreo del conjunto principal ------------------------------------
    print("\n[ muestreando conjunto principal... ]")
    secs, longs, neff_real, descart = sampling.muestrear_banda(
        BANDA_PRINCIPAL, cfg["n_eff_principal"], rng_py, min_L=LAGS[-1] + 1
    )
    v_concat = np.concatenate(secs)
    linea(f"- trayectorias válidas: {len(secs):,} (descartadas por L<{LAGS[-1]+1}: {descart})")
    linea(f"- N_eff(lag1) alcanzado: {neff_real:,}")
    linea(f"- L: media={longs.mean():.1f}, mediana={np.median(longs):.0f}, "
          f"min={longs.min()}, max={longs.max()}")

    # -- A. Sanity check ----------------------------------------------------
    linea("\n## A. Sanity check")
    ev = float(v_concat.mean())
    drift, _ = stats.drift_geometrico(cfg["n_drift"], rng=np.random.default_rng(cfg["semilla"] + 2))
    linea(f"- E[v] (sobre {len(v_concat):,} pasos del conjunto principal) = {ev:.5f}  (teoría 2.0)")
    linea(f"- drift geométrico (muestreo independiente)               = {drift:.5f}  (teoría 0.75)")
    linea("- distribución de v vs 2^-k:")
    linea("  | v | empírico | 2^-k |")
    linea("  |---|----------|------|")
    for k in range(1, 7):
        emp = float((v_concat == k).mean())
        linea(f"  | {k} | {emp:.4f} | {2.0**-k:.4f} |")

    # -- B/C. ACF agregado bruto + corregido --------------------------------
    print("\n[ ACF agregado y corrección de sesgo... ]")
    agg = acf.acf_agregado(secs, max_lag=LAGS[-1])
    linea("\n## B / C. ACF bruto, sesgo -1/(L-1) y ACF corregido")
    linea("| lag | bruto | baseline | corregido | banda ±1.96/√Neff | z | N_eff |")
    linea("|-----|-------|----------|-----------|-------------------|---|-------|")
    for k in LAGS:
        linea(f"| {k} | {fmt(agg['rho_bruto'][k])} | {fmt(agg['baseline_medio'][k])} | "
              f"{fmt(agg['rho_corregido'][k])} | ±{agg['banda'][k]:.4f} | "
              f"{agg['z'][k]:+.2f} | {int(agg['n_eff'][k]):,} |")
    linea("\nInterpretación: el sesgo de muestra finita empuja TODOS los lags hacia "
          "negativo (~-1/(L-1)). El ACF corregido descuenta ese efecto; |z|>1.96 es "
          "el criterio analítico (provisional: supone colas livianas).")

    # -- D. Permutación (principal) -----------------------------------------
    print("\n[ prueba de permutación (principal)... esto tarda ]")
    perm = acf.prueba_permutacion(secs, LAGS, B=cfg["B"], rng=rng_np)
    linea("\n## D. Prueba de PERMUTACIÓN (principal)")
    linea("| lag | observado | nulo medio | nulo [2.5%, 97.5%] | p(1 cola sup) | p(2 colas) | signif. |")
    linea("|-----|-----------|------------|--------------------|---------------|------------|---------|")
    for k in LAGS:
        d = perm[k]
        linea(f"| {k} | {fmt(d['observado'])} | {fmt(d['nulo_media'])} | "
              f"[{d['p2_5']:+.4f}, {d['p97_5']:+.4f}] | {d['p_una_cola_sup']:.4f} | "
              f"{d['p_dos_colas']:.4f} | {'sí' if d['significativo'] else 'no'} |")
    linea("\nLa permutación NO asume colas livianas: es la prueba que manda. El nulo "
          "ya incorpora el sesgo -1/(L-1) y la cola pesada de v.")

    veredicto_analitico = abs(agg["z"][1]) > 1.96 and agg["rho_corregido"][1] > 0
    veredicto_perm = perm[1]["significativo"] and perm[1]["observado"] > perm[1]["nulo_media"]
    if veredicto_analitico != veredicto_perm:
        linea(f"\n⚠ Discrepancia analítica vs permutación en lag-1 "
              f"(analítica={'sí' if veredicto_analitico else 'no'}, "
              f"permutación={'sí' if veredicto_perm else 'no'}). Prevalece la permutación.")

    # -- E. Estratificación por L -------------------------------------------
    print("\n[ estratificación por longitud L... ]")
    linea("\n## E. Estratificación por longitud L (confound de supervivencia)")
    linea("| bucket | n | L medio | lag1 bruto | lag1 corregido | perm lag1 (p 1 cola) | signif. |")
    linea("|--------|---|---------|------------|----------------|----------------------|---------|")
    grupos = sampling.estratificar_por_longitud(secs)
    lag1_buckets = []  # (etiqueta, corregido, half_width, significativo, n)
    for etiqueta, _, _ in sampling.BUCKETS_L:
        g = grupos[etiqueta]
        if len(g) == 0:
            linea(f"| {etiqueta} | 0 | — | — | — | — | — |")
            continue
        agg_g = acf.acf_agregado(g, max_lag=1)
        lmed = float(np.mean([len(x) for x in g]))
        corr = agg_g["rho_corregido"][1]
        if len(g) >= MIN_TRAYECTORIAS_BUCKET:
            perm_g = acf.prueba_permutacion(g, [1], B=cfg["B"], rng=rng_np)[1]
            half = (perm_g["p97_5"] - perm_g["p2_5"]) / 2.0
            sig = perm_g["significativo"] and perm_g["observado"] > perm_g["nulo_media"]
            pstr = f"{perm_g['p_una_cola_sup']:.4f}"
            sigstr = "sí" if sig else "no"
            lag1_buckets.append((etiqueta, corr, half, sig, len(g)))
        else:
            pstr, sigstr = "n<min", "—"
        linea(f"| {etiqueta} | {len(g):,} | {lmed:.1f} | {fmt(agg_g['rho_bruto'][1])} | "
              f"{fmt(corr)} | {pstr} | {sigstr} |")
    linea("\nObjetivo: si el lag-1 solo aparece al mezclar largos pero muere DENTRO de "
          "cada estrato, sería paradoja de Simpson (artefacto de supervivencia).")

    # -- F. Escalamiento por bits -------------------------------------------
    print("\n[ escalamiento por bits (motor de enteros Python)... ]")
    linea("\n## F. Escalamiento por bits (¿lag-1 estable o se desvanece?)")
    linea("| banda bits | n tray | L medio | N_eff | lag1 bruto | lag1 corregido | perm (p 1 cola) | signif. |")
    linea("|------------|--------|---------|-------|------------|----------------|-----------------|---------|")
    lag1_bits = []  # (etiqueta, corregido, half, sig)
    for banda in BANDAS_BITS:
        s_b, l_b, neff_b, _ = sampling.muestrear_banda(
            banda, cfg["n_eff_bits"], rng_py, min_L=LAGS[-1] + 1
        )
        agg_b = acf.acf_agregado(s_b, max_lag=1)
        perm_b = acf.prueba_permutacion(s_b, [1], B=cfg["B"], rng=rng_np)[1]
        corr = agg_b["rho_corregido"][1]
        half = (perm_b["p97_5"] - perm_b["p2_5"]) / 2.0
        sig = perm_b["significativo"] and perm_b["observado"] > perm_b["nulo_media"]
        etiqueta = f"{banda[0]}-{banda[1]}"
        lag1_bits.append((etiqueta, corr, half, sig))
        linea(f"| {etiqueta} | {len(s_b):,} | {l_b.mean():.1f} | {neff_b:,} | "
              f"{fmt(agg_b['rho_bruto'][1])} | {fmt(corr)} | {perm_b['p_una_cola_sup']:.4f} | "
              f"{'sí' if sig else 'no'} |")
    linea("\nObjetivo: si el lag-1 corregido se mantiene ~estable al subir bits, apunta "
          "a estructura 2-ádica real; si se desvanece hacia 0, era artefacto de números chicos.")

    # -- Figuras ------------------------------------------------------------
    print("\n[ guardando figuras... ]")
    figs = []
    figs.append(viz.acf_con_banda(agg, "rho_bruto", "outputs/it2_acf_bruto.png",
                                  "ACF bruto vs lag (banda ±1.96/√Neff)"))
    figs.append(viz.acf_con_banda(agg, "rho_corregido", "outputs/it2_acf_corregido.png",
                                  "ACF corregido (−1/(L−1) descontado) vs lag"))
    figs.append(viz.histograma_longitudes(longs))
    figs.append(viz.distribucion_v(v_concat))
    if lag1_buckets:
        figs.append(viz.lag1_por_grupo(
            [b[0] for b in lag1_buckets], [b[1] for b in lag1_buckets],
            [b[2] for b in lag1_buckets], [b[2] for b in lag1_buckets],
            "outputs/it2_lag1_por_L.png", "Lag-1 corregido por bucket de L",
            "bucket de longitud L"))
    figs.append(viz.lag1_por_grupo(
        [b[0] for b in lag1_bits], [b[1] for b in lag1_bits],
        [b[2] for b in lag1_bits], [b[2] for b in lag1_bits],
        "outputs/it2_lag1_por_bits.png", "Lag-1 corregido por banda de bits",
        "banda de bits"))
    linea("\n## Figuras")
    for f in figs:
        linea(f"- {f}")

    # -- G. Veredicto -------------------------------------------------------
    c1 = bool(agg["rho_corregido"][1] > 0)
    c2 = bool(perm[1]["significativo"] and perm[1]["observado"] > perm[1]["nulo_media"])
    buckets_validos = [b for b in lag1_buckets]
    c3 = bool(buckets_validos) and all(b[1] > 0 and b[3] for b in buckets_validos)
    todos_pos_sig = all(b[1] > 0 and b[3] for b in lag1_bits)
    val_lo = next((b[1] for b in lag1_bits if b[0] == "20-24"), None)
    val_hi = next((b[1] for b in lag1_bits if b[0] == "60-64"), None)
    no_colapsa = (val_lo is not None and val_hi is not None and val_lo > 0
                  and val_hi >= 0.5 * val_lo)
    c4 = bool(todos_pos_sig and no_colapsa)

    linea("\n## G. VEREDICTO")
    linea(f"- (1) lag-1 corregido positivo:                 {'SÍ' if c1 else 'NO'} "
          f"({fmt(agg['rho_corregido'][1])})")
    linea(f"- (2) significativo por permutación (positivo): {'SÍ' if c2 else 'NO'} "
          f"(p 1 cola = {perm[1]['p_una_cola_sup']:.4f})")
    linea(f"- (3) estable DENTRO de buckets de L:           {'SÍ' if c3 else 'NO'}")
    linea(f"- (4) estable al escalar bits (no colapsa):     {'SÍ' if c4 else 'NO'}")
    lags_negativos = [k for k in LAGS[1:] if perm[k]["significativo"]
                      and perm[k]["observado"] < perm[k]["nulo_media"]]
    linea(f"- lags 2..10 significativos negativos tras permutación: "
          f"{lags_negativos if lags_negativos else 'ninguno (compatibles con sesgo)'}")

    sobrevive = c1 and c2 and c3 and c4
    linea("\n### Regla de decisión")
    if sobrevive:
        linea("**Las CUATRO condiciones se cumplen.** El lag-1 positivo sobrevive a la "
              "corrección de sesgo, a la permutación, a la estratificación por longitud y "
              "al escalamiento por bits. Bajo la regla acordada, la iteración 3 queda "
              "justificada. Lenguaje prudente: es señal empírica robusta de dependencia de "
              "corto alcance, NO una afirmación matemática probada.")
    else:
        fallidas = [n for n, c in zip("1234", [c1, c2, c3, c4]) if not c]
        linea(f"**Falla(n) la(s) condición(es) {', '.join(fallidas)}.** Bajo la regla "
              "acordada, NO se justifica la iteración 3 por esta vía: el ACF de lag-1 "
              "aparente es atribuible a artefacto de muestra finita / supervivencia / "
              "números chicos. Se documenta y se CIERRA la veta sin romantizarla.")

    linea("\n### Pregunta de interpretación para la iteración 3 (anotada, NO resuelta)")
    linea("Si el lag-1 resultara real: ¿es reducible a la mecánica determinista del mapa "
          "más la distribución marginal de v (y por tanto 'esperable'), o hay algo no "
          "trivial encima? No se formula nada sobre Tao todavía.")
    linea("\n" + "=" * 70)

    # -- Escribir reporte ---------------------------------------------------
    docs = Path(__file__).resolve().parents[1] / "docs"
    docs.mkdir(exist_ok=True)
    reporte = docs / "iteracion2_reporte.md"
    encabezado = ["# Iteración 2 — ¿sobrevive el lag-1 positivo?", "",
                  f"_Generado por `scripts/run_iteration2.py` (semilla {cfg['semilla']}"
                  f"{', modo rápido' if args.rapido else ''})._", ""]
    reporte.write_text("\n".join(encabezado + R), encoding="utf-8")
    print(f"\nReporte escrito en: {reporte}")


if __name__ == "__main__":
    main()
