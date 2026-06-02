# Experimento confirmatorio pre-cierre — bits vs L

_Generado por `scripts/run_confirmatorio.py` (semilla 20260615)._

======================================================================
  Collatz — experimento confirmatorio pre-cierre (NO iteración 3)
  Des-confundir magnitud (bits) de largo (L) + asíntota del lag-1
======================================================================

## Regla de decisión (FIJADA antes de calcular)
Dos vías independientes de des-confusión: (i) lag-1 a L FIJA por banda de bits; (ii) regresión parcial del lag-1 sobre (bits, L), coeficiente de bits.
- **RATIFICAR CIERRE** si, a L fija, el lag-1 sigue decayendo con la magnitud: coef. de bits con IC95 enteramente < 0 **y** el ajuste de asíntota favorece '→ 0' (o la constante c tiene IC que incluye 0). Documentar como efecto de tamaño finito / magnitud y cerrar con honor.
- **NO CERRAR (abrir iteración 3 SOLO reducibilidad)** si, a L fija, el lag-1 queda plano: coef. de bits NO significativo, **o** asíntota c claramente positiva (IC95 enteramente > 0). No perseguir magnitud, no tocar Tao. Registrar que el escenario más probable sigue siendo estructura 2-ádica conocida.

## Parámetros
- semilla: 20260615 | B permutaciones: 1000
- bandas de bits: ['20-24', '30-34', '40-44', '50-54', '60-64', '70-74', '80-84']
- trayectorias por banda: 20,000
- umbral N_eff para ventana de L común: 25,000
- enteros nativos de Python (sin overflow)

## Muestreo por banda
| banda | n tray | L medio | L min | L max |
|-------|--------|---------|-------|-------|
| 20-24 | 20,000 | 53.6 | 11 | 166 |
| 30-34 | 20,000 | 77.8 | 11 | 282 |
| 40-44 | 20,000 | 101.6 | 16 | 338 |
| 50-54 | 20,000 | 125.6 | 29 | 341 |
| 60-64 | 20,000 | 149.8 | 38 | 354 |
| 70-74 | 20,000 | 173.9 | 51 | 416 |
| 80-84 | 20,000 | 198.0 | 71 | 457 |

## 3. Des-confusión por ventana de L FIJA
Ventana de L común más ancha hallada: **[86, 136)** (umbral N_eff = 25,000).

| banda | n en ventana | L medio | N_eff | lag1 corregido | banda ±1.96/√Neff | perm p(2 colas) | signif. |
|-------|--------------|---------|-------|----------------|-------------------|-----------------|---------|
| 20-24 | 1,946 | 98.8 | 190,222 | +0.0511 | ±0.0045 | 0.0010 | sí |
| 30-34 | 6,570 | 103.7 | 674,818 | +0.0364 | ±0.0024 | 0.0010 | sí |
| 40-44 | 10,220 | 107.7 | 1,090,344 | +0.0263 | ±0.0019 | 0.0010 | sí |
| 50-54 | 10,132 | 111.3 | 1,117,359 | +0.0181 | ±0.0019 | 0.0010 | sí |
| 60-64 | 7,103 | 114.8 | 808,518 | +0.0108 | ±0.0022 | 0.0010 | sí |
| 70-74 | 3,789 | 118.3 | 444,435 | +0.0097 | ±0.0029 | 0.0010 | sí |
| 80-84 | 1,623 | 121.0 | 194,774 | +0.0063 | ±0.0044 | 0.0030 | sí |

Dispersión de L media dentro de la ventana: 22.3 pasos (min 98.8, max 121.0). Si es chica frente a la L típica, L está bien controlada y la comparación entre bandas es ~apples-to-apples; el confound residual lo arbitra la regresión parcial (3bis), que es el des-confound robusto.

## 3bis. Des-confusión por regresión parcial (todos los datos)
- n trayectorias en la regresión: 140,000
- coef. **bits** (efecto de magnitud A IGUAL L): -0.000624 por bit (IC95 [-0.000659, -0.000589])
- coef. **L** (efecto del largo a igual magnitud):  +0.000067 por paso (IC95 [+0.000056, +0.000077])
- coef. de bits significativamente negativo: sí
Interpretación: el signo y significancia de 'bits' dice si la magnitud sigue moviendo el lag-1 una vez fijado L. 'L' dice lo análogo para el largo.

## 4. Ajuste de asíntota del lag-1
Serie ajustada: lag-1 a L FIJA.
- modelo (a) → 0:        RSS=5.38e-06  AIC=-94.55
- modelo (b) → c:         RSS=5.47e-06  AIC=-92.43
  asíntota c = -0.00191  (IC95 [-0.01684, +0.01301])
- modelo favorecido por AIC: **cero**

## 5. Estructura multi-lag que sobrevive a la permutación (banda 40-44)
- trayectorias usadas: 5,000
| lag | observado | nulo medio | corregido | signif. |
|-----|-----------|------------|-----------|---------|
| 1 | +0.0120 | -0.0098 | +0.0219 | sí |
| 2 | -0.0243 | -0.0097 | -0.0147 | sí |
| 3 | -0.0002 | -0.0095 | +0.0093 | sí |
| 4 | +0.0316 | -0.0094 | +0.0411 | sí |
| 5 | +0.0082 | -0.0093 | +0.0174 | sí |
| 6 | -0.0071 | -0.0092 | +0.0021 | no |
| 7 | +0.0121 | -0.0091 | +0.0213 | sí |
| 8 | +0.0042 | -0.0090 | +0.0132 | sí |
| 9 | -0.0228 | -0.0089 | -0.0140 | sí |
| 10 | -0.0146 | -0.0087 | -0.0059 | sí |
- lags significativos: 9/10
- fracción de energía ACF que sobrevive a la permutación: 0.682 (p=0.0010)
Lectura prudente (sin interpretar de más): la secuencia de v's es determinista dada n, así que una fracción alta que sobrevive a la permutación intra-trayectoria es estructura dependiente del ORDEN, es decir, mecánica del mapa.

## Figuras
- outputs/conf_lag1_fijoL.png
- outputs/conf_multilag.png

## 6. VEREDICTO (según la regla pre-registrada)
- ¿coef. de bits (a L fija) significativamente negativo? SÍ
- ¿asíntota constante claramente positiva (IC95>0)?      NO
- ¿lag-1 a L fija plano (coef. bits no significativo)?    NO

**RATIFICAR CIERRE.** A L fija el lag-1 sigue decayendo con la magnitud y la asíntota es compatible con 0. El lag-1 aparente es un efecto de tamaño finito / magnitud, no estructura persistente. Se cierra la veta del lag-1 con honor. Lenguaje prudente: esto es evidencia empírica de un artefacto de tamaño, no una demostración.

======================================================================