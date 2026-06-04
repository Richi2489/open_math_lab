"""
run_explicit_formula_demo.py — Exhibición B: los ceros reconstruyen los primos
==============================================================================

Visualización pedagógica de la fórmula explícita de Riemann (1859): cómo la suma sobre
los ceros de ζ corrige la densidad suave de primos hasta dibujar la escalera ψ(x). NO es
un intento de RH; es la "segunda cara" del objeto (la §3 fue la textura estadística GUE,
esto es el contenido aritmético).

Escribe docs/explicit_formula.md y figuras en docs/figures/riemann_explicit/.
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

from riemann import explicit_formula as ef  # noqa: E402
from riemann import viz, zeros  # noqa: E402

FIGDIR = "docs/figures/riemann_explicit"
NS = [10, 50, 200, 1000]
X_MIN, X_MAX = 2.0, 100.0


def main():
    R = []

    def linea(s=""):
        print(s)
        R.append(s)

    linea("=" * 72)
    linea("  Riemann Lab — Exhibición B: los ceros reconstruyen los primos")
    linea("  (fórmula explícita de Riemann, 1859; NO es un intento de RH)")
    linea("=" * 72)

    gammas = zeros.cargar_ceros(max(NS))
    x = np.linspace(X_MIN, X_MAX, 3000)
    psi_ex = ef.psi_exacta(x, X_MAX)
    curvas = [(N, ef.psi_aprox(x, gammas[:N])) for N in NS]

    # -- 1. La idea --------------------------------------------------------
    linea("\n## 1. La idea")
    linea("ψ(x) = Σ_{p^k ≤ x} ln(p) cuenta (con peso) las potencias de primo: una escalera "
          "que salta ln(p) en cada p^k. La fórmula explícita la reconstruye como")
    linea("    ψ_N(x) = x − 2√x·Σ_{n≤N} Re( x^{iγ_n}/(½+iγ_n) ) − ln(2π) − ½ln(1−x⁻²).")
    linea("El término `x` es la tendencia suave (teorema de los números primos); cada cero "
          "γ_n aporta una FRECUENCIA que corrige esa tendencia hacia los saltos reales.")

    # -- 2. Convergencia ---------------------------------------------------
    linea("\n## 2. Convergencia con el número de ceros")
    linea("| N ceros | error absoluto medio en [2,100] |")
    linea("|---------|----------------------------------|")
    for N, ap in curvas:
        linea(f"| {N} | {float(np.mean(np.abs(ap - psi_ex))):.4f} |")
    linea("\nCon pocos ceros es una onda burda; al crecer N las ondas se enganchan a los "
          "saltos en los primos y el error medio decrece.")

    # -- 3. Sanity PNT -----------------------------------------------------
    r = ef.psi_exacta(1e5) / 1e5
    linea("\n## 3. Sanity check (teorema de los números primos)")
    linea(f"- ψ(10⁵)/10⁵ = {r:.5f}  (→ 1, como debe).")

    # -- 4. Gibbs ----------------------------------------------------------
    xz = np.linspace(6.0, 14.0, 1500)
    psi_z = ef.psi_exacta(xz, X_MAX)
    curvas_z = [(N, ef.psi_aprox(xz, gammas[:N])) for N in NS]
    over = float(np.max(curvas_z[-1][1] - psi_z))
    linea("\n## 4. Honestidad: el fenómeno de Gibbs")
    linea(f"- Cerca de cada salto hay overshoot/oscilación; con N={NS[-1]} el sobre-impulso "
          f"máximo en [6,14] es ≈ {over:.2f}.")
    linea("- Es el comportamiento ESPERADO de una suma truncada de Fourier (fenómeno de "
          "Gibbs), NO un error del método ni evidencia de nada sobre RH.")

    # -- Figuras -----------------------------------------------------------
    figs = [
        viz.formula_explicita_panel(x, psi_ex, curvas, f"{FIGDIR}/psi_reconstruccion.png"),
        viz.formula_explicita_gibbs(xz, psi_z, curvas_z, f"{FIGDIR}/psi_gibbs.png"),
    ]
    linea("\n## Figuras")
    for fg in figs:
        linea(f"- {fg}")

    # -- Cierre ------------------------------------------------------------
    linea("\n## La segunda cara del objeto")
    linea("La §3 del lab (estadística de espaciados, GUE) mostró la TEXTURA del espectro de "
          "ζ. Esto muestra su CONTENIDO ARITMÉTICO: los mismos ceros, leídos como "
          "frecuencias, reconstruyen el conteo de primos. Densidad suave, estadística GUE y "
          "huella de los primos son tres caras del mismo objeto — y por eso un operador que "
          "las unifique sería tan codiciado. Lenguaje prudente: esto es la visualización de "
          "un resultado clásico (la fórmula explícita), no evidencia sobre RH.")
    linea("\n" + "=" * 72)

    out = Path(__file__).resolve().parents[1] / "docs" / "explicit_formula.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    cab = ["# La fórmula explícita — los ceros reconstruyen los primos", "",
           "_Demostración expositiva del Riemann Statistical Lab. NO es un intento de RH._",
           f"_Generado por `scripts/run_explicit_formula_demo.py` (N∈{NS}, x∈[{X_MIN:.0f},"
           f"{X_MAX:.0f}])._", ""]
    out.write_text("\n".join(cab + R), encoding="utf-8")
    print(f"\nReporte escrito en: {out}")


if __name__ == "__main__":
    main()
