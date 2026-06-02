"""
run_confirmatorio.py — des-confundir magnitud (bits) de largo (L) antes de cerrar
=================================================================================

Experimento confirmatorio pre-cierre (NO es iteración 3). La iteración 2 cerró la
veta del lag-1 porque decaía con los bits, pero ese test confunde magnitud con largo.
Aquí se des-confunde y se caracteriza la asíntota, con la regla de decisión FIJADA
de antemano (se imprime antes de calcular nada).

Salidas: docs/confirmatorio_reporte.md y figuras en outputs/.

Uso:
    python scripts/run_confirmatorio.py            # completo
    python scripts/run_confirmatorio.py --rapido   # chico, para fumar el pipeline
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

from collatz import acf, deconfound, sampling, viz  # noqa: E402

BANDAS_BITS = [(20, 24), (30, 34), (40, 44), (50, 54), (60, 64), (70, 74), (80, 84)]
LAGS = list(range(1, 11))


def etiqueta(banda):
    return f"{banda[0]}-{banda[1]}"


def centro(banda):
    return (banda[0] + banda[1]) / 2.0


def fmt(x, n=4):
    return f"{x:+.{n}f}"


def construir_config(rapido, semilla):
    if rapido:
        return dict(semilla=semilla, B=300, n_por_banda=2500,
                    umbral_neff=5000, n_multilag=1000)
    return dict(semilla=semilla, B=1000, n_por_banda=20000,
                umbral_neff=25000, n_multilag=5000)


def main():
    parser = argparse.ArgumentParser(description="Confirmatorio: des-confundir bits vs L.")
    parser.add_argument("--rapido", action="store_true")
    parser.add_argument("--semilla", type=int, default=20260615)
    args = parser.parse_args()
    cfg = construir_config(args.rapido, args.semilla)

    rng_py = random.Random(cfg["semilla"])
    rng_np = np.random.default_rng(cfg["semilla"] + 1)
    R = []

    def linea(s=""):
        print(s)
        R.append(s)

    linea("=" * 70)
    linea("  Collatz — experimento confirmatorio pre-cierre (NO iteración 3)")
    linea("  Des-confundir magnitud (bits) de largo (L) + asíntota del lag-1")
    linea("=" * 70)

    # -- Regla de decisión PRE-REGISTRADA (antes de ver números) ------------
    linea("\n## Regla de decisión (FIJADA antes de calcular)")
    linea("Dos vías independientes de des-confusión: (i) lag-1 a L FIJA por banda de "
          "bits; (ii) regresión parcial del lag-1 sobre (bits, L), coeficiente de bits.")
    linea("- **RATIFICAR CIERRE** si, a L fija, el lag-1 sigue decayendo con la magnitud: "
          "coef. de bits con IC95 enteramente < 0 **y** el ajuste de asíntota favorece "
          "'→ 0' (o la constante c tiene IC que incluye 0). Documentar como efecto de "
          "tamaño finito / magnitud y cerrar con honor.")
    linea("- **NO CERRAR (abrir iteración 3 SOLO reducibilidad)** si, a L fija, el lag-1 "
          "queda plano: coef. de bits NO significativo, **o** asíntota c claramente "
          "positiva (IC95 enteramente > 0). No perseguir magnitud, no tocar Tao. Registrar "
          "que el escenario más probable sigue siendo estructura 2-ádica conocida.")

    # -- Parámetros ---------------------------------------------------------
    linea("\n## Parámetros")
    linea(f"- semilla: {cfg['semilla']} | B permutaciones: {cfg['B']}")
    linea(f"- bandas de bits: {[etiqueta(b) for b in BANDAS_BITS]}")
    linea(f"- trayectorias por banda: {cfg['n_por_banda']:,}")
    linea(f"- umbral N_eff para ventana de L común: {cfg['umbral_neff']:,}")
    linea(f"- enteros nativos de Python (sin overflow)")

    # -- Muestreo de las 7 bandas (guardando secuencias) --------------------
    print("\n[ muestreando 7 bandas con enteros de Python... ]")
    seqs_por_banda = {}
    centros = {}
    linea("\n## Muestreo por banda")
    linea("| banda | n tray | L medio | L min | L max |")
    linea("|-------|--------|---------|-------|-------|")
    for b in BANDAS_BITS:
        seqs, longs, _ = sampling.muestrear_banda_por_conteo(
            b, cfg["n_por_banda"], rng_py, min_L=LAGS[-1] + 1
        )
        seqs_por_banda[etiqueta(b)] = seqs
        centros[etiqueta(b)] = centro(b)
        linea(f"| {etiqueta(b)} | {len(seqs):,} | {longs.mean():.1f} | "
              f"{longs.min()} | {longs.max()} |")

    # -- Búsqueda de ventana de L común -------------------------------------
    print("\n[ buscando ventana de L común a todas las bandas... ]")
    ventana = deconfound.buscar_ventana_L(seqs_por_banda, cfg["umbral_neff"])
    linea("\n## 3. Des-confusión por ventana de L FIJA")
    if ventana is None:
        linea(f"⚠ No existe ventana de L con N_eff>={cfg['umbral_neff']:,} en TODAS las "
              "bandas (magnitud y L están demasiado correlacionadas en los extremos). "
              "Se omite el análisis a L fija; se confía en la regresión parcial (sección 3bis).")
        fijo = None
    else:
        lo, hi = ventana
        linea(f"Ventana de L común más ancha hallada: **[{lo}, {hi})** "
              f"(umbral N_eff = {cfg['umbral_neff']:,}).")
        fijo = deconfound.lag1_por_banda_en_ventana(seqs_por_banda, lo, hi, cfg["B"], rng_np)
        linea("\n| banda | n en ventana | L medio | N_eff | lag1 corregido | banda ±1.96/√Neff | perm p(2 colas) | signif. |")
        linea("|-------|--------------|---------|-------|----------------|-------------------|-----------------|---------|")
        for b in BANDAS_BITS:
            e = etiqueta(b)
            d = fijo[e]
            linea(f"| {e} | {d['n']:,} | {d['L_medio']:.1f} | {d['n_eff']:,} | "
                  f"{fmt(d['lag1_corregido'])} | ±{d['banda_signif']:.4f} | "
                  f"{d['perm']['p_dos_colas']:.4f} | "
                  f"{'sí' if d['perm']['significativo'] else 'no'} |")
        l_medios = [fijo[etiqueta(b)]["L_medio"] for b in BANDAS_BITS]
        spread = max(l_medios) - min(l_medios)
        linea(f"\nDispersión de L media dentro de la ventana: {spread:.1f} pasos "
              f"(min {min(l_medios):.1f}, max {max(l_medios):.1f}). "
              "Si es chica frente a la L típica, L está bien controlada y la comparación "
              "entre bandas es ~apples-to-apples; el confound residual lo arbitra la "
              "regresión parcial (3bis), que es el des-confound robusto.")

    # -- Regresión parcial (de-confound robusto) ----------------------------
    print("\n[ regresión parcial lag-1 ~ bits + L ... ]")
    reg = deconfound.regresion_parcial(seqs_por_banda, centros)
    b_bits = reg["beta"][1]
    ic_bits = reg["ic95"][1]
    b_L = reg["beta"][2]
    ic_L = reg["ic95"][2]
    linea("\n## 3bis. Des-confusión por regresión parcial (todos los datos)")
    linea(f"- n trayectorias en la regresión: {reg['n']:,}")
    linea(f"- coef. **bits** (efecto de magnitud A IGUAL L): {b_bits:+.6f} por bit "
          f"(IC95 [{ic_bits[0]:+.6f}, {ic_bits[1]:+.6f}])")
    linea(f"- coef. **L** (efecto del largo a igual magnitud):  {b_L:+.6f} por paso "
          f"(IC95 [{ic_L[0]:+.6f}, {ic_L[1]:+.6f}])")
    bits_sig_neg = ic_bits[1] < 0
    bits_no_sig = ic_bits[0] <= 0 <= ic_bits[1]
    linea(f"- coef. de bits significativamente negativo: {'sí' if bits_sig_neg else 'no'}")
    linea("Interpretación: el signo y significancia de 'bits' dice si la magnitud sigue "
          "moviendo el lag-1 una vez fijado L. 'L' dice lo análogo para el largo.")

    # -- Ajuste de asíntota (sobre el lag-1 a L fija si existe) -------------
    print("\n[ ajustando asíntota... ]")
    linea("\n## 4. Ajuste de asíntota del lag-1")
    if fijo is not None:
        xs = [centros[etiqueta(b)] for b in BANDAS_BITS]
        ys = [fijo[etiqueta(b)]["lag1_corregido"] for b in BANDAS_BITS]
        sig = [max(fijo[etiqueta(b)]["banda_signif"], 1e-6) for b in BANDAS_BITS]
        fuente = "lag-1 a L FIJA"
    else:
        # fallback: usar el lag-1 crudo por banda (confundido, pero algo es algo)
        xs, ys, sig = [], [], []
        for b in BANDAS_BITS:
            agg = acf.acf_agregado(seqs_por_banda[etiqueta(b)], max_lag=1)
            xs.append(centro(b)); ys.append(float(agg["rho_corregido"][1]))
            sig.append(max(float(agg["banda"][1]), 1e-6))
        fuente = "lag-1 por banda (sin L fija — confundido)"
    ajuste = deconfound.ajustar_asintota(xs, ys, sigma=sig)
    linea(f"Serie ajustada: {fuente}.")
    mc, mk = ajuste.get("modelo_cero", {}), ajuste.get("modelo_const", {})
    if "aic" in mc:
        linea(f"- modelo (a) → 0:        RSS={mc.get('rss', float('nan')):.2e}  AIC={mc['aic']:.2f}")
    if "aic" in mk:
        c = mk.get("asintota", float("nan")); se = mk.get("se_asintota", float("nan"))
        ic_lo, ic_hi = c - 1.96 * se, c + 1.96 * se
        linea(f"- modelo (b) → c:         RSS={mk.get('rss', float('nan')):.2e}  AIC={mk['aic']:.2f}")
        linea(f"  asíntota c = {c:+.5f}  (IC95 [{ic_lo:+.5f}, {ic_hi:+.5f}])")
    linea(f"- modelo favorecido por AIC: **{ajuste['ganador']}**")

    asintota_pos = False
    if "asintota" in mk and np.isfinite(mk.get("se_asintota", np.nan)):
        c = mk["asintota"]; se = mk["se_asintota"]
        asintota_pos = (ajuste["ganador"] == "const") and (c - 1.96 * se > 0)

    # -- Supervivencia multi-lag --------------------------------------------
    print("\n[ supervivencia multi-lag a la permutación... ]")
    secs_ml = seqs_por_banda["40-44"][: cfg["n_multilag"]]
    ml = deconfound.supervivencia_multilag(secs_ml, LAGS, cfg["B"], rng_np)
    linea("\n## 5. Estructura multi-lag que sobrevive a la permutación (banda 40-44)")
    linea(f"- trayectorias usadas: {len(secs_ml):,}")
    linea("| lag | observado | nulo medio | corregido | signif. |")
    linea("|-----|-----------|------------|-----------|---------|")
    for k in LAGS:
        d = ml["por_lag"][k]
        linea(f"| {k} | {fmt(d['observado'])} | {fmt(d['nulo_media'])} | "
              f"{fmt(d['corregido'])} | {'sí' if d['significativo'] else 'no'} |")
    linea(f"- lags significativos: {ml['lags_significativos']}/{ml['n_lags']}")
    linea(f"- fracción de energía ACF que sobrevive a la permutación: "
          f"{ml['fraccion_energia_sobrevive']:.3f} (p={ml['p_energia']:.4f})")
    linea("Lectura prudente (sin interpretar de más): la secuencia de v's es determinista "
          "dada n, así que una fracción alta que sobrevive a la permutación intra-trayectoria "
          "es estructura dependiente del ORDEN, es decir, mecánica del mapa.")

    # -- Figuras ------------------------------------------------------------
    print("\n[ guardando figuras... ]")
    figs = []
    figs.append(viz.lag1_fijo_L_con_ajustes(
        xs, ys, sig, ajuste, "outputs/conf_lag1_fijoL.png",
        f"lag-1(bits) — {fuente}"))
    figs.append(viz.supervivencia_multilag(ml["por_lag"]))
    linea("\n## Figuras")
    for f in figs:
        linea(f"- {f}")

    # -- 6. Veredicto pre-registrado ----------------------------------------
    ratificar = bool(bits_sig_neg and (ajuste["ganador"] == "cero" or not asintota_pos))
    linea("\n## 6. VEREDICTO (según la regla pre-registrada)")
    linea(f"- ¿coef. de bits (a L fija) significativamente negativo? {'SÍ' if bits_sig_neg else 'NO'}")
    linea(f"- ¿asíntota constante claramente positiva (IC95>0)?      {'SÍ' if asintota_pos else 'NO'}")
    linea(f"- ¿lag-1 a L fija plano (coef. bits no significativo)?    {'SÍ' if bits_no_sig else 'NO'}")
    if ratificar:
        linea("\n**RATIFICAR CIERRE.** A L fija el lag-1 sigue decayendo con la magnitud y "
              "la asíntota es compatible con 0. El lag-1 aparente es un efecto de tamaño "
              "finito / magnitud, no estructura persistente. Se cierra la veta del lag-1 con "
              "honor. Lenguaje prudente: esto es evidencia empírica de un artefacto de "
              "tamaño, no una demostración.")
    else:
        linea("\n**NO CERRAR LIMPIAMENTE.** A L fija el lag-1 NO decae claramente a 0 (o la "
              "asíntota es positiva). El decaimiento por bits de la iteración 2 estaba "
              "confundido con L. Bajo la regla pre-registrada NO se persigue la magnitud ni "
              "se toca Tao: si acaso, la iteración 3 se abre SOLO sobre la pregunta de "
              "reducibilidad (¿el patrón multi-lag completo se explica por mecánica del mapa "
              "+ marginal de v?). Escenario más probable: estructura 2-ádica conocida, no "
              "señal nueva.")
    linea("\n" + "=" * 70)

    docs = Path(__file__).resolve().parents[1] / "docs"
    docs.mkdir(exist_ok=True)
    reporte = docs / "confirmatorio_reporte.md"
    encabezado = ["# Experimento confirmatorio pre-cierre — bits vs L", "",
                  f"_Generado por `scripts/run_confirmatorio.py` (semilla {cfg['semilla']}"
                  f"{', modo rápido' if args.rapido else ''})._", ""]
    reporte.write_text("\n".join(encabezado + R), encoding="utf-8")
    print(f"\nReporte escrito en: {reporte}")


if __name__ == "__main__":
    main()
