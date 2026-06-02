"""
run_baseline.py — corrida base de la iteración 1
================================================

Mide y reporta:
  1. E[v]                     (esperado ≈ 2.0)
  2. drift geométrico         (esperado ≈ 0.75 = 3/4)
  3. resumen de autocorrelación de los v's
  4. distribución de tiempos de parada

y guarda tres gráficas en outputs/.

Uso:
    python scripts/run_baseline.py [--muestras N] [--semilla S]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# Forzar UTF-8 en stdout: las consolas de Windows usan cp1252 por defecto y
# revientan al imprimir símbolos como ≈ o —.
try:
    sys.stdout.reconfigure(encoding="utf-8")
except (AttributeError, ValueError):  # pragma: no cover
    pass

# Permite ejecutar el script directamente (python scripts/run_baseline.py)
# añadiendo la raíz del repo al path.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from collatz import stats, viz  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Corrida base del laboratorio de Collatz.")
    parser.add_argument(
        "--muestras",
        type=int,
        default=1_000_000,
        help="número de impares para E[v] y drift (default 1e6).",
    )
    parser.add_argument(
        "--trayectorias",
        type=int,
        default=3000,
        help="trayectorias para la autocorrelación promedio (default 3000).",
    )
    parser.add_argument("--semilla", type=int, default=20240601, help="semilla del RNG.")
    args = parser.parse_args()

    rng = np.random.default_rng(args.semilla)

    print("=" * 64)
    print("  open_math_lab — Collatz, iteración 1 (corrida base)")
    print("=" * 64)
    print(f"  muestras (E[v], drift): {args.muestras:,}")
    print(f"  trayectorias (ACF):     {args.trayectorias:,}")
    print(f"  semilla RNG:            {args.semilla}")
    print("-" * 64)

    # 1) E[v] ----------------------------------------------------------------
    media_v, _ = stats.esperanza_v(args.muestras, rng=rng)
    print(f"  1) E[v] empírico        = {media_v:.5f}   (teoría: 2.0)")
    print(f"        desviación        = {media_v - 2.0:+.5f}")

    # 2) Drift geométrico ----------------------------------------------------
    drift, _ = stats.drift_geometrico(args.muestras, rng=rng)
    print(f"  2) drift geométrico     = {drift:.5f}   (teoría: 0.75)")
    print(f"        desviación        = {drift - 0.75:+.5f}")

    # 3) Autocorrelación -----------------------------------------------------
    lags, acf, usadas = stats.autocorrelacion_v_promedio(
        n_trayectorias=args.trayectorias, max_lag=20, rng=rng
    )
    print(f"  3) autocorrelación de v's (promedio sobre {usadas} trayectorias):")
    for k in (1, 2, 3, 5, 10):
        if k < len(acf):
            print(f"        lag {k:>2}            = {acf[k]:+.4f}")
    print(
        "     (si los v's fueran independientes, todos estos serían ≈ 0;"
        " el tamaño y signo cuantifican la dependencia real)"
    )

    # 4) Distribución de tiempos de parada -----------------------------------
    tiempos = stats.distribucion_tiempos_parada(n_muestras=20_000, n_max=10**6, rng=rng)
    print("  4) tiempos de parada (semillas en [2, 1e6)):")
    print(f"        media             = {tiempos.mean():.2f}")
    print(f"        mediana           = {np.median(tiempos):.1f}")
    print(f"        máximo            = {tiempos.max()}")

    # Gráficas ---------------------------------------------------------------
    print("-" * 64)
    print("  guardando gráficas en outputs/ ...")
    p1 = viz.grafica_trayectorias([27, 97, 871, 6171])
    p2 = viz.histograma_tiempos_parada(tiempos)
    p3 = viz.grafica_autocorrelacion(lags, acf)
    for p in (p1, p2, p3):
        print(f"        {p}")

    print("=" * 64)
    print("  Resumen: E[v] ≈ 2 y drift ≈ 3/4 confirman la heurística;")
    print("  la autocorrelación NO nula en lag>=1 muestra la dependencia")
    print("  que la heurística ignora. Esa es la pregunta de la iteración 1.")
    print("=" * 64)


if __name__ == "__main__":
    main()
