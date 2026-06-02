# Confirmación alta-T con ceros de Odlyzko — iteración 2

_Generado por `scripts/run_odlyzko_highT.py` (N=2000000, unfold=local_density, GUE 2000×2000×40, semilla 2026)._

========================================================================
  Riemann Lab — confirmación alta-T con ceros de Odlyzko (iteración 2)
========================================================================

## A. Regla pre-registrada (con la trampa del KS)
Comparar SIEMPRE contra GUE de tamaño finito simulado, no contra la fórmula asintótica. **Trampa del KS:** con ~2M de muestras el KS es hipersensible y casi seguro dará p<0.05; eso NO es anomalía, es la corrección conocida de altura finita al GUE. El veredicto se lee por la TENDENCIA (¿la desviación encoge vs baja altura y por ventanas?), por la separación clara frente a Poisson, y por si el residuo cabe en la variabilidad del GUE finito — NUNCA por el p-valor del KS.

## B. Dataset y validaciones
- archivo: `data/odlyzko/zeros6` | total en archivo: 2,001,052 ceros
- muestra usada: N = 2,000,000 (start=0)
- altura T: 14.135 … 1131944.5
- validación: 1D numérico, estrictamente ascendente, sin duplicados/NaN/inf ✓

## C. Método de unfolding
- usado: **local_density**  (s_n = (γ_{n+1}−γ_n)·log(γ_n/2π)/2π)

## D. Sanity checks
- media del espaciado unfolded (debe ≈ 1): 1.00000
- cross-check Odlyzko vs mpmath (primeros 2000 γ): máx |Δ| = 5.00e-10 (precisión Odlyzko ~4e-9) ✓
- GUE finito (2000×2000, 63,960 espac.): var=0.1819, P(s<0.5)=0.1149

## E. ζ alta-T vs GUE, Poisson, Wigner
- ζ (N=1,999,999): media=1.0000, var=0.1661, P(s<0.5)=0.1021, P(s<0.25)=0.0139
| comparación | KS | KS p-valor | Wasserstein |
|-------------|----|-----------|-------------|
| ζ vs GUE | 0.0165 | 4.47e-15 | 0.0163 |
| ζ vs Poisson | 0.2926 | 0.00e+00 | 0.4289 |
- KS p-valor ζ-GUE = 4.47e-15: con N=1,999,999 esto es la TRAMPA esperada (hipersensibilidad), no una anomalía. Lo informativo es el TAMAÑO de la desviación (KS=0.0165, Wasserstein=0.0163) y su tendencia.

## F. Baja altura vs Odlyzko alta-T
| dataset | N spac | T_min | T_max | var | P(s<0.5) | KS vs GUE | Wass vs GUE | KS vs Poisson |
|---------|--------|-------|-------|-----|----------|-----------|-------------|---------------|
| mpmath N=2000 | 1,999 | 14 | 2515 | 0.1477 | 0.0800 | 0.0442 | 0.0349 | 0.3158 |
| Odlyzko primeros 2000 | 1,999 | 14 | 2515 | 0.1477 | 0.0800 | 0.0442 | 0.0349 | 0.3158 |
| Odlyzko N=2,000,000 | 1,999,999 | 14 | 1131944 | 0.1661 | 0.1021 | 0.0165 | 0.0163 | 0.2926 |

(GUE finito de referencia: var=0.1819, P(s<0.5)=0.1149. Leer cómo |var−var_GUE| y |P(s<0.5)−P_GUE| cambian al subir la altura.)

## G. Estabilidad por ventanas no solapadas (altura creciente)
| ventana | T_min | T_max | N spac | var | P(s<0.5) | |Δvar| | KS vs GUE | Wass vs GUE |
|---------|-------|-------|--------|-----|----------|--------|-----------|-------------|
| W1 | 14 | 260877 | 399,999 | 0.1635 | 0.0993 | 0.0184 | 0.0195 | 0.0191 |
| W2 | 260877 | 489737 | 399,999 | 0.1659 | 0.1019 | 0.0160 | 0.0167 | 0.0166 |
| W3 | 489738 | 709043 | 399,999 | 0.1666 | 0.1024 | 0.0153 | 0.0160 | 0.0158 |
| W4 | 709043 | 922554 | 399,999 | 0.1670 | 0.1030 | 0.0148 | 0.0152 | 0.0153 |
| W5 | 922555 | 1131944 | 399,999 | 0.1674 | 0.1038 | 0.0145 | 0.0154 | 0.0149 |

- tendencia |Δvar| con altura: pendiente -4.02e-09 (W1=0.0184 → W5=0.0145)
- tendencia Wasserstein con altura: pendiente -4.38e-09 (W1=0.0191 → W5=0.0149)

## Figuras
- docs/figures/riemann_odlyzko/hist_spacing.png
- docs/figures/riemann_odlyzko/cdf.png
- docs/figures/riemann_odlyzko/cdf_diferencia.png
- docs/figures/riemann_odlyzko/por_ventana.png
- docs/figures/riemann_odlyzko/baja_vs_alta.png
- docs/figures/riemann_odlyzko/r2.png

## H. Limitaciones
- Altura ~10⁶, no ~10²⁰: los datasets de muy alta altura (offsets) usan otro formato/parser; quedan para una iteración futura.
- GUE finito de dim fija como referencia (cuasi-asintótica para el spacing NN); no se modela la correspondencia exacta dim↔altura.
- KS hipersensible a N grande (trampa pre-registrada); por eso el veredicto se apoya en tamaño de desviación y tendencia, no en el p-valor.
- R₂ con unfolding local (cumsum de espaciados), suficiente para la firma de repulsión, no para correlaciones de largo alcance finas.

## I. Veredicto (lenguaje prudente)
- ¿la desviación a GUE ENCOGE vs baja altura? SÍ (|Δvar|: 0.0342 a alturas ~10³ → 0.0158 a ~10⁶)
- ¿Poisson claramente descartado? SÍ (KS ζ-Poisson=0.2926 ≫ KS ζ-GUE=0.0165)
- ¿GUE sigue siendo el baseline correcto? SÍ (var y P(s<0.5) de ζ ≈ GUE finito)
- ¿desviación estable/decreciente por ventanas (no residuo creciente)? SÍ
- KS p-valor ζ-GUE = 4.5e-15: trampa esperada, NO se interpreta como anomalía.

**SE CONFIRMA EL CIERRE DE LA RÉPLICA.** A altura ~10⁶ los ceros de ζ están mucho más cerca del GUE finito que de Poisson, y la pequeña desviación de baja altura ENCOGE al subir la altura (consistente con la corrección conocida de tamaño/altura finita). No hay evidencia de anomalía estructural; el p-valor del KS es la trampa de muestra grande, no señal. La repulsión tipo GUE de los gaps unfolded de ζ queda reproducida limpiamente. Prudencia: es una réplica empírica de un resultado conocido (Montgomery–Odlyzko), no evidencia sobre RH.

========================================================================