# Iteración 2 — ¿sobrevive el lag-1 positivo?

_Generado por `scripts/run_iteration2.py` (semilla 20260601)._

======================================================================
  open_math_lab — Collatz, iteración 2
  ¿sobrevive el lag-1 positivo tras corrección Y permutación?
======================================================================

## 1. Reproducibilidad / parámetros
- semilla maestra: 20260601
- permutaciones B: 1000
- N_eff objetivo (banda principal (40, 44) bits): 300,000
- N_eff objetivo por banda de bits: 200,000
- lags analizados: 1..10
- min_L (para soportar todos los lags): 11
- enteros: nativos de Python (sin riesgo de overflow)
- trayectorias válidas: 2,969 (descartadas por L<11: 0)
- N_eff(lag1) alcanzado: 300,034
- L: media=102.1, mediana=99, min=24, max=311

## A. Sanity check
- E[v] (sobre 303,003 pasos del conjunto principal) = 1.99422  (teoría 2.0)
- drift geométrico (muestreo independiente)               = 0.74852  (teoría 0.75)
- distribución de v vs 2^-k:
  | v | empírico | 2^-k |
  |---|----------|------|
  | 1 | 0.4993 | 0.5000 |
  | 2 | 0.2464 | 0.2500 |
  | 3 | 0.1244 | 0.1250 |
  | 4 | 0.0727 | 0.0625 |
  | 5 | 0.0308 | 0.0312 |
  | 6 | 0.0132 | 0.0156 |

## B / C. ACF bruto, sesgo -1/(L-1) y ACF corregido
| lag | bruto | baseline | corregido | banda ±1.96/√Neff | z | N_eff |
|-----|-------|----------|-----------|-------------------|---|-------|
| 1 | +0.0154 | -0.0099 | +0.0253 | ±0.0036 | +13.87 | 300,034 |
| 2 | -0.0243 | -0.0099 | -0.0144 | ±0.0036 | -7.86 | 297,065 |
| 3 | +0.0001 | -0.0099 | +0.0100 | ±0.0036 | +5.40 | 294,096 |
| 4 | +0.0311 | -0.0099 | +0.0409 | ±0.0036 | +22.09 | 291,127 |
| 5 | +0.0073 | -0.0098 | +0.0172 | ±0.0037 | +9.21 | 288,158 |
| 6 | -0.0092 | -0.0098 | +0.0006 | ±0.0037 | +0.33 | 285,189 |
| 7 | +0.0149 | -0.0098 | +0.0247 | ±0.0037 | +13.14 | 282,220 |
| 8 | +0.0108 | -0.0098 | +0.0206 | ±0.0037 | +10.89 | 279,251 |
| 9 | -0.0251 | -0.0098 | -0.0153 | ±0.0037 | -8.07 | 276,282 |
| 10 | -0.0178 | -0.0098 | -0.0080 | ±0.0037 | -4.20 | 273,313 |

Interpretación: el sesgo de muestra finita empuja TODOS los lags hacia negativo (~-1/(L-1)). El ACF corregido descuenta ese efecto; |z|>1.96 es el criterio analítico (provisional: supone colas livianas).

## D. Prueba de PERMUTACIÓN (principal)
| lag | observado | nulo medio | nulo [2.5%, 97.5%] | p(1 cola sup) | p(2 colas) | signif. |
|-----|-----------|------------|--------------------|---------------|------------|---------|
| 1 | +0.0154 | -0.0098 | [-0.0132, -0.0066] | 0.0010 | 0.0010 | sí |
| 2 | -0.0243 | -0.0097 | [-0.0132, -0.0064] | 1.0000 | 0.0010 | sí |
| 3 | +0.0001 | -0.0095 | [-0.0128, -0.0064] | 0.0010 | 0.0010 | sí |
| 4 | +0.0311 | -0.0094 | [-0.0129, -0.0061] | 0.0010 | 0.0010 | sí |
| 5 | +0.0073 | -0.0093 | [-0.0126, -0.0060] | 0.0010 | 0.0010 | sí |
| 6 | -0.0092 | -0.0093 | [-0.0128, -0.0060] | 0.4985 | 0.9690 | no |
| 7 | +0.0149 | -0.0089 | [-0.0120, -0.0058] | 0.0010 | 0.0010 | sí |
| 8 | +0.0108 | -0.0088 | [-0.0121, -0.0056] | 0.0010 | 0.0010 | sí |
| 9 | -0.0251 | -0.0088 | [-0.0120, -0.0055] | 1.0000 | 0.0010 | sí |
| 10 | -0.0178 | -0.0087 | [-0.0121, -0.0056] | 1.0000 | 0.0010 | sí |

La permutación NO asume colas livianas: es la prueba que manda. El nulo ya incorpora el sesgo -1/(L-1) y la cola pesada de v.

## E. Estratificación por longitud L (confound de supervivencia)
| bucket | n | L medio | lag1 bruto | lag1 corregido | perm lag1 (p 1 cola) | signif. |
|--------|---|---------|------------|----------------|----------------------|---------|
| L<25 | 1 | 24.0 | -0.0279 | +0.0156 | n<min | — |
| 25-50 | 94 | 41.4 | -0.0190 | +0.0058 | n<min | — |
| 50-100 | 1,404 | 78.3 | +0.0032 | +0.0161 | 0.0010 | sí |
| 100-200 | 1,453 | 127.5 | +0.0231 | +0.0310 | 0.0010 | sí |
| L>=200 | 17 | 224.6 | +0.0289 | +0.0333 | n<min | — |

Objetivo: si el lag-1 solo aparece al mezclar largos pero muere DENTRO de cada estrato, sería paradoja de Simpson (artefacto de supervivencia).

## F. Escalamiento por bits (¿lag-1 estable o se desvanece?)
| banda bits | n tray | L medio | N_eff | lag1 bruto | lag1 corregido | perm (p 1 cola) | signif. |
|------------|--------|---------|-------|------------|----------------|-----------------|---------|
| 20-24 | 3,778 | 53.9 | 200,006 | +0.0329 | +0.0518 | 0.0010 | sí |
| 30-34 | 2,633 | 77.0 | 200,029 | +0.0189 | +0.0320 | 0.0010 | sí |
| 40-44 | 2,013 | 100.4 | 200,015 | +0.0089 | +0.0190 | 0.0010 | sí |
| 50-54 | 1,579 | 127.7 | 200,097 | +0.0084 | +0.0163 | 0.0010 | sí |
| 60-64 | 1,349 | 149.3 | 200,066 | +0.0059 | +0.0127 | 0.0010 | sí |

Objetivo: si el lag-1 corregido se mantiene ~estable al subir bits, apunta a estructura 2-ádica real; si se desvanece hacia 0, era artefacto de números chicos.

## Figuras
- outputs/it2_acf_bruto.png
- outputs/it2_acf_corregido.png
- outputs/it2_longitudes.png
- outputs/it2_dist_v.png
- outputs/it2_lag1_por_L.png
- outputs/it2_lag1_por_bits.png

## G. VEREDICTO
- (1) lag-1 corregido positivo:                 SÍ (+0.0253)
- (2) significativo por permutación (positivo): SÍ (p 1 cola = 0.0010)
- (3) estable DENTRO de buckets de L:           SÍ
- (4) estable al escalar bits (no colapsa):     NO
- lags 2..10 significativos negativos tras permutación: [2, 9, 10]

### Regla de decisión
**Falla(n) la(s) condición(es) 4.** Bajo la regla acordada, NO se justifica la iteración 3 por esta vía: el ACF de lag-1 aparente es atribuible a artefacto de muestra finita / supervivencia / números chicos. Se documenta y se CIERRA la veta sin romantizarla.

### Pregunta de interpretación para la iteración 3 (anotada, NO resuelta)
Si el lag-1 resultara real: ¿es reducible a la mecánica determinista del mapa más la distribución marginal de v (y por tanto 'esperable'), o hay algo no trivial encima? No se formula nada sobre Tao todavía.

======================================================================