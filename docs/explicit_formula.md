# La fórmula explícita — los ceros reconstruyen los primos

_Demostración expositiva del Riemann Statistical Lab. NO es un intento de RH._
_Generado por `scripts/run_explicit_formula_demo.py` (N∈[10, 50, 200, 1000], x∈[2,100])._

========================================================================
  Riemann Lab — Exhibición B: los ceros reconstruyen los primos
  (fórmula explícita de Riemann, 1859; NO es un intento de RH)
========================================================================

## 1. La idea
ψ(x) = Σ_{p^k ≤ x} ln(p) cuenta (con peso) las potencias de primo: una escalera que salta ln(p) en cada p^k. La fórmula explícita la reconstruye como
    ψ_N(x) = x − 2√x·Σ_{n≤N} Re( x^{iγ_n}/(½+iγ_n) ) − ln(2π) − ½ln(1−x⁻²).
El término `x` es la tendencia suave (teorema de los números primos); cada cero γ_n aporta una FRECUENCIA que corrige esa tendencia hacia los saltos reales.

## 2. Convergencia con el número de ceros
| N ceros | error absoluto medio en [2,100] |
|---------|----------------------------------|
| 10 | 0.7532 |
| 50 | 0.4121 |
| 200 | 0.1963 |
| 1000 | 0.0707 |

Con pocos ceros es una onda burda; al crecer N las ondas se enganchan a los saltos en los primos y el error medio decrece.

## 3. Sanity check (teorema de los números primos)
- ψ(10⁵)/10⁵ = 1.00052  (→ 1, como debe).

## 4. Honestidad: el fenómeno de Gibbs
- Cerca de cada salto hay overshoot/oscilación; con N=1000 el sobre-impulso máximo en [6,14] es ≈ 0.98.
- Es el comportamiento ESPERADO de una suma truncada de Fourier (fenómeno de Gibbs), NO un error del método ni evidencia de nada sobre RH.

## Figuras
- docs/figures/riemann_explicit/psi_reconstruccion.png
- docs/figures/riemann_explicit/psi_gibbs.png

## La segunda cara del objeto
La §3 del lab (estadística de espaciados, GUE) mostró la TEXTURA del espectro de ζ. Esto muestra su CONTENIDO ARITMÉTICO: los mismos ceros, leídos como frecuencias, reconstruyen el conteo de primos. Densidad suave, estadística GUE y huella de los primos son tres caras del mismo objeto — y por eso un operador que las unifique sería tan codiciado. Lenguaje prudente: esto es la visualización de un resultado clásico (la fórmula explícita), no evidencia sobre RH.

========================================================================