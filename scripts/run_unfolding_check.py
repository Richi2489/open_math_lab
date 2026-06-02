"""
run_unfolding_check.py — descargar la regla pre-registrada (barato, sin descargas)
==================================================================================

La iteración 1 reprodujo la repulsión tipo GUE, pero ζ vs GUE no era exactamente cero
(KS ≈ 0.042, p ≈ 0.01): a baja altura ζ reprime un pelín más que el GUE finito. La regla
anti-autoengaño exige comprobar que eso NO se debe al unfolding y que se comporta como un
efecto de altura finita, ANTES de interpretarlo. Aquí se hace con los ~2000 ceros ya
cacheados (sin Odlyzko todavía):

  1. UNFOLDING EXACTO: re-unfold con θ(T)/π+1 (Riemann–Siegel) en vez de la asintótica;
     re-medir KS / Wasserstein / P(s<0.5). ¿Cambia algo? (Se espera sub-permille.)
  2. TENDENCIA POR ALTURA: cuartiles de γ; ¿la brecha ζ−GUE encoge al subir la altura
     DENTRO de nuestra propia muestra?

Veredicto por la regla pre-registrada. Lenguaje prudente.
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

from riemann import gue, metrics, spacing, viz, zeros  # noqa: E402

N_CEROS = 2000
SEMILLA = 20260602
DIM_GUE = 1500
N_MATRICES = 30  # ~36k espaciados de referencia


def f(x, n=4):
    return f"{x:.{n}f}"


def main():
    rng = np.random.default_rng(SEMILLA)
    R = []

    def linea(s=""):
        print(s)
        R.append(s)

    linea("=" * 70)
    linea("  Riemann Lab — descarga empírica de la regla pre-registrada")
    linea("=" * 70)
    linea("\n## Regla (recordatorio)")
    linea("Una desviación de GUE NO se interpreta si desaparece al cambiar unfolding, "
          "ventana de altura, tamaño de muestra o baseline. Aquí se prueban unfolding y "
          "altura con la muestra ya cacheada (sin descargas).")

    print("\n[ cargando ceros + unfolding exacto (siegeltheta, cacheado)... ]")
    gammas = zeros.cargar_ceros(N_CEROS)
    zeros.validar_ceros(gammas)
    w_asint = spacing.unfold(gammas)
    w_exact = zeros.cargar_unfold_exacto(N_CEROS, verbose=True)
    s_asint = np.diff(w_asint)
    s_exact = np.diff(w_exact)

    # Referencia GUE finito grande (cuasi-asintótica) y Poisson para contexto.
    s_gue = gue.espaciados_gue(DIM_GUE, rng, n_matrices=N_MATRICES, recorte=0.1)
    var_gue = float(s_gue.var())
    p_gue = float((s_gue < 0.5).mean())

    # ------------------------------------------------------------------
    # 1. Unfolding exacto vs asintótico
    # ------------------------------------------------------------------
    linea("\n## 1. Unfolding exacto (θ/π+1) vs asintótico (N_suave)")
    linea(f"- diferencia máxima en la coordenada unfolded |w_exact − w_asint| = "
          f"{np.max(np.abs(w_exact - w_asint)):.2e}")
    linea(f"- diferencia máxima en los espaciados |s_exact − s_asint| = "
          f"{np.max(np.abs(s_exact - s_asint)):.2e}")
    linea("\n| métrica (vs GUE) | asintótico | exacto | cambio |")
    linea("|------------------|-----------|--------|--------|")
    filas = []
    for nombre, s in [("asint", s_asint), ("exacto", s_exact)]:
        d = metrics.distancias(s, s_gue)
        m = metrics.momentos(s)
        filas.append((d["ks"], d["wasserstein"], m["frac_menor_0.5"], m["var"]))
    (ks_a, wa_a, p_a, v_a), (ks_e, wa_e, p_e, v_e) = filas
    linea(f"| KS | {f(ks_a)} | {f(ks_e)} | {ks_e - ks_a:+.2e} |")
    linea(f"| Wasserstein | {f(wa_a)} | {f(wa_e)} | {wa_e - wa_a:+.2e} |")
    linea(f"| P(s<0.5) | {f(p_a)} | {f(p_e)} | {p_e - p_a:+.2e} |")
    linea(f"| var | {f(v_a)} | {f(v_e)} | {v_e - v_a:+.2e} |")
    unfold_estable = abs(ks_e - ks_a) < 5e-3 and abs(v_e - v_a) < 1e-3
    linea(f"\n→ unfolding como causa: {'DESCARTADO (cambio ínfimo)' if unfold_estable else 'NO descartable'}.")

    # ------------------------------------------------------------------
    # 2. Tendencia por cuartiles de altura (con unfolding EXACTO)
    # ------------------------------------------------------------------
    linea("\n## 2. Tendencia por altura — cuartiles de γ (unfolding exacto)")
    linea("| cuartil | rango T | n | var(ζ) | P(s<0.5) ζ | |Δvar| | KS vs GUE |")
    linea("|---------|---------|---|--------|------------|--------|-----------|")
    grupos = np.array_split(np.arange(N_CEROS), 4)
    centros_T, var_gaps, ks_list = [], [], []
    for q, idx in enumerate(grupos, 1):
        wq = w_exact[idx]
        sq = np.diff(wq)
        vq = float(sq.var())
        pq = float((sq < 0.5).mean())
        dvar = abs(vq - var_gue)
        ksq = float(metrics.distancias(sq, s_gue)["ks"])
        Tmin, Tmax = float(gammas[idx[0]]), float(gammas[idx[-1]])
        centros_T.append(0.5 * (Tmin + Tmax))
        var_gaps.append(dvar)
        ks_list.append(ksq)
        linea(f"| Q{q} | {Tmin:.0f}–{Tmax:.0f} | {len(sq):,} | {f(vq)} | {f(pq)} | "
              f"{f(dvar)} | {f(ksq)} |")
    se_var = var_gue * np.sqrt(2.0 / len(grupos[0]))
    linea(f"\n- referencia GUE: var = {var_gue:.4f}, P(s<0.5) = {p_gue:.4f} "
          f"(matriz {DIM_GUE}×{DIM_GUE}, {len(s_gue):,} espaciados).")
    linea(f"- error estándar aprox. de var(ζ) por cuartil ≈ ±{se_var:.4f} "
          "(potencia limitada con ~500 espaciados/cuartil; se lee la TENDENCIA, no cada punto).")

    # Tendencia: ¿encoge la brecha de varianza con la altura?
    encoge_var = var_gaps[-1] < var_gaps[0]
    pend_var = float(np.polyfit(centros_T, var_gaps, 1)[0])
    pend_ks = float(np.polyfit(centros_T, ks_list, 1)[0])
    linea(f"\n- brecha |Δvar|: Q1 = {var_gaps[0]:.4f} → Q4 = {var_gaps[-1]:.4f} "
          f"(pendiente {pend_var:+.2e}/unidad T)")
    linea(f"- KS(ζ vs GUE): Q1 = {ks_list[0]:.4f} → Q4 = {ks_list[-1]:.4f} "
          f"(pendiente {pend_ks:+.2e}/unidad T)")

    # ------------------------------------------------------------------
    # Figura
    # ------------------------------------------------------------------
    fig = viz.tendencia_altura(centros_T, var_gaps, ks_list)
    linea(f"\n## Figura\n- {fig}")

    # ------------------------------------------------------------------
    # Veredicto
    # ------------------------------------------------------------------
    linea("\n## Veredicto (regla pre-registrada)")
    linea(f"- estable bajo unfolding exacto: {'SÍ' if unfold_estable else 'NO'}")
    linea(f"- brecha de varianza encoge con la altura (Q4<Q1): {'SÍ' if encoge_var else 'NO'} "
          f"(pendiente var {pend_var:+.1e}, KS {pend_ks:+.1e})")
    if unfold_estable and encoge_var:
        linea("\n**REGLA DESCARGADA.** La pequeña desviación ζ−GUE NO se debe al unfolding "
              "(cambio ínfimo con θ exacta) y ENCOGE con la altura dentro de la propia "
              "muestra: se comporta como un efecto de tamaño finito / baja altura, no como "
              "una desviación estructural. No amerita teoría para 'explicarla'. Prudencia: "
              "esto es consistencia empírica con altura finita, no una demostración.")
    elif unfold_estable and not encoge_var:
        linea("\n**PARCIAL.** El unfolding queda descartado como causa, pero la brecha NO "
              "encoge claramente con la altura en esta muestra (potencia limitada a ~2000 "
              "ceros). Se ANOTA para confirmación futura con tablas de Odlyzko de alta T, "
              "donde la altura finita se separa limpiamente. No se interpreta aún.")
    else:
        linea("\n**REVISAR.** El resultado depende del unfolding — antes de cualquier "
              "lectura hay que estabilizar el método de unfolding.")
    linea("\n" + "=" * 70)

    docs = Path(__file__).resolve().parents[1] / "docs"
    docs.mkdir(exist_ok=True)
    rep = docs / "unfolding_height_check_report.md"
    cab = ["# Descarga de la regla pre-registrada — unfolding exacto + tendencia por altura",
           "", f"_Generado por `scripts/run_unfolding_check.py` "
           f"(N={N_CEROS}, semilla {SEMILLA})._", ""]
    rep.write_text("\n".join(cab + R), encoding="utf-8")
    print(f"\nReporte escrito en: {rep}")


if __name__ == "__main__":
    main()
