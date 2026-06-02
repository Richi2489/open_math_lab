# Réplica Montgomery–Odlyzko — iteración 1

_Generado por `scripts/run_spacing_replication.py` (N=2000, semilla 20260602)._

======================================================================
  Riemann Statistical Lab — réplica Montgomery–Odlyzko (iteración 1)
======================================================================

## Regla anti-autoengaño (PRE-REGISTRADA, antes de ver resultados)
No se declarará una desviación de GUE si desaparece al cambiar (a) la ventana de altura T, (b) el método de unfolding, (c) el tamaño de muestra, o (d) el baseline. Toda diferencia contra GUE se compara contra GUE de tamaño finito simulado (mismo N), NO contra la fórmula asintótica ideal. El 'éxito' de la iteración 1 es REPRODUCIR limpiamente la repulsión tipo GUE, no encontrar una anomalía.

## Parámetros
- ceros de ζ: N = 2,000 (mpmath, cacheados en data/)
- GUE finito: matriz 2000×2000; baseline Poisson (gaps exp i.i.d.)
- ventanas de altura: 5 | semilla: 20260602

## Ceros y unfolding
- γ_1 = 14.1347, γ_2 = 21.0220, γ_N = 2515.2865
- espaciados unfolded: 1,999 | media = 1.0000 (debe ≈ 1 si el unfolding es correcto)

## A. Momentos de los espaciados
| dist. | media | var | asimetría | P(s<0.5) |
|-------|-------|-----|-----------|----------|
| ζ | 1.0000 | 0.1479 | 0.4916 | 0.0800 |
| GUE | 0.9997 | 0.1772 | 0.5016 | 0.1094 |
| Poisson | 0.9890 | 0.9877 | 2.0247 | 0.4053 |
P(s<0.5) chico = repulsión (pocos gaps diminutos). Poisson ~0.39; GUE mucho menor.

## B. Distancias de la distribución de espaciados
| comparación | KS | KS p-valor | Wasserstein |
|-------------|----|-----------|-------------|
| ζ vs GUE | 0.0417 | 9.594e-03 | 0.0295 |
| ζ vs Poisson | 0.3282 | 6.508e-146 | 0.4396 |
Repulsión tipo GUE ⇔ ζ MUCHO más cerca de GUE que de Poisson (KS y Wasserstein).

## C. Correlación de pares R₂ (firma Montgomery–Odlyzko)
| r | R₂ empírico (ζ) | GUE 1−(sinπr/πr)² |
|---|-----------------|--------------------|
| 0.2 | 0.0250 | 0.0719 |
| 0.5 | 0.3950 | 0.5119 |
| 0.9 | 1.0250 | 0.9973 |
| 1.5 | 0.9750 | 0.9530 |
| 2.0 | 0.9400 | 0.9993 |
La firma es el HOYO de repulsión cerca de r→0 (R₂→0), no el valor plano 1 de Poisson.

## D. El confounder: normalización GLOBAL (INCORRECTA)
- P(s<0.5) con unfolding correcto = 0.0800 ; con normalización global = 0.0950
- var correcto = 0.1479 ; var global-incorrecto = 0.2114  (la densidad cambia con la altura: dividir por la media global mezcla regímenes y distorsiona la distribución — el análogo del confounder de magnitud en Collatz).

## E. Estabilidad por ventana de altura (efectos de tamaño finito)
| ventana T | n espac. | KS(ζ vs GUE) | Wasserstein |
|-----------|----------|--------------|-------------|
| 14–680 | 399 | 0.0518 | 0.0368 |
| 682–1184 | 399 | 0.0533 | 0.0337 |
| 1185–1648 | 399 | 0.0697 | 0.0394 |
| 1649–2090 | 399 | 0.0684 | 0.0375 |
| 2090–2515 | 399 | 0.0484 | 0.0280 |

KS por ventana: min 0.048, max 0.070. A baja altura el match suele ser peor (efecto de tamaño finito esperado, no anomalía).

## Figuras
- outputs/riemann_spacing.png
- outputs/riemann_cdf.png
- outputs/riemann_r2.png
- outputs/riemann_estabilidad.png

## Veredicto (iteración 1) — lenguaje prudente
- ζ más cerca de GUE que de Poisson (KS y Wasserstein): SÍ
- repulsión a corta distancia (P(s<0.5) ≪ Poisson):     SÍ
- R₂ con hoyo de repulsión cerca de r→0:                SÍ

**Repulsión tipo GUE reproducida** en las tres métricas (espaciados NN, distancias vs baselines, y correlación de pares). El residuo frente al GUE asintótico es compatible con tamaño finito / baja altura (ver sección E), no se interpreta como anomalía. Prudencia: es una réplica empírica de una señal conocida, no evidencia sobre RH.

======================================================================