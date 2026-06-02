# ADR 0001 — Cierre de la veta del lag-1 (dinámica impar de Collatz)

- **Estado:** Aceptado
- **Fecha:** 2026-06-01
- **Commit que congela el estado experimental:** `f57c9abd50ca2eca36d3e2da47acfbcab6989dfc`
- **Decisión:** RATIFICAR EL CIERRE de la hipótesis del lag-1 como señal nueva.

---

## Contexto

El laboratorio mide empíricamente la heurística probabilística de Collatz sobre el mapa
acelerado de impares `T(n) = (3n+1)/2^v`. La iteración 1 confirmó las dos predicciones
del modelo (E[v] ≈ 2.0, drift geométrico ≈ 3/4) y observó una autocorrelación de **lag-1
positiva** (≈ +0.0137 en bruto) en la secuencia de los `v`.

**Hipótesis evaluada:** ese lag-1 positivo es una correlación de corto alcance *real* en
la dinámica, es decir, estructura no capturada por el modelo de `v` i.i.d.

La heurística asume independencia entre pasos; cualquier dependencia robusta y persistente
sería el punto de partida interesante. Por eso la hipótesis merecía una prueba falsable
seria antes de aceptarse o descartarse.

## Regla de decisión (pre-registrada, fijada ANTES de ver los números del confirmatorio)

El lag-1 se acepta como señal estructural persistente **solo si** sobrevive a las cuatro:

1. ser positivo tras corregir el sesgo de muestra finita `-1/(L-1)`;
2. ser significativo bajo la prueba de **permutación** intra-trayectoria (B ≥ 1000);
3. mantenerse dentro de los estratos de longitud `L` (no ser paradoja de Simpson);
4. mantenerse estable al escalar la **magnitud** (bits del arranque), no desvanecerse.

Confirmatorio adicional (des-confound), también pre-registrado:

- **RATIFICAR CIERRE** si, a **L fija**, el lag-1 sigue decayendo con la magnitud
  (coef. de bits con IC95 enteramente < 0) **y** el ajuste de asíntota favorece "→ 0"
  (o la constante `c` tiene IC95 que incluye 0).
- **NO CERRAR** (abrir iteración 3 solo sobre reducibilidad) si a L fija el lag-1 queda
  plano (coef. de bits no significativo) **o** la asíntota `c` es claramente positiva.

## Resultados

**Iteración 2** (semilla 20260601, B=1000):
- lag-1 corregido = **+0.0253** (z ≈ +13.9), permutación p = 0.0010 → pasa (1), (2), (3).
- Escalamiento por bits: lag-1 corregido decae monótonamente **+0.052 → +0.013**
  (20-24 → 60-64 bits) → **falla (4)**. Veredicto provisional: cerrar. *Pero* ese test
  confunde magnitud con longitud `L`.

**Confirmatorio** (semilla 20260615, B=1000, 20k trayectorias/banda):
- Ventana de L fija y angosta **[86, 136)** (L media 99–121 entre bandas, dispersión 22):
  a mismo largo, el lag-1 corregido **aún decae con la magnitud** (+0.0511 → +0.0063).
- Regresión parcial `lag1 ~ bits + L` (n = 140,000): coef. de **bits = −0.000624/bit**
  (IC95 [−0.000659, −0.000589]) controlando L; coef. de **L = +0.000067/paso** (IC95 > 0).
  Dos efectos reales y opuestos; a igual L la magnitud empuja el lag-1 hacia 0.
- Asíntota: modelo "→ 0" gana por AIC (−94.55 vs −92.43); constante alternativa
  c = −0.0019 (IC95 [−0.0168, +0.0130], incluye 0).

## Decisión

**RATIFICAR EL CIERRE.** Las dos vías de des-confound coinciden: a longitud fija el lag-1
sigue decayendo con la magnitud y la asíntota es estadísticamente compatible con cero. El
lag-1 observado es un **efecto de tamaño finito**, no evidencia de correlación estructural
persistente.

## Consecuencias

- La veta del lag-1 queda **cerrada y documentada** (post-mortem en
  `docs/collatz_lag1_closure.md`). Resultado negativo, reproducible.
- **No** se persigue una demostración de reducibilidad del patrón multi-lag residual: solo
  confirmaría estructura 2-ádica ya conocida (vector de paridad, Terras 1976) y no aportaría
  descubrimiento. Queda anotado, no abierto.
- No se formula nada sobre la conjetura ni sobre líneas de ataque tipo Tao. Lenguaje
  prudente: esto es evidencia empírica, no una demostración.
- El valor del ejercicio es **metodológico**: heurística → hipótesis → falsación con regla
  pre-registrada → des-confound (incluido un confounder detectado y corregido a mitad de
  camino). Reutilizable en otros sistemas.
