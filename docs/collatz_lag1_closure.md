# Post-mortem — la veta del lag-1 en la dinámica impar de Collatz

> Documento de cierre. Cuenta el arco completo y honesto de una hipótesis que se
> investigó con rigor y se descartó. Es un **resultado negativo bien hecho**: el valor
> está en el método, no en un descubrimiento. Lenguaje prudente a propósito: nada aquí
> es una demostración, ni toca la conjetura ni líneas de ataque tipo Tao.

---

## 1. El punto de partida: la heurística calza

El mapa acelerado sobre impares es `T(n) = (3n+1) / 2^v`, con `v = v₂(3n+1) ≥ 1`. La
heurística probabilística modela cada `v` como geométrica i.i.d. con `P(v=k) = 2^-k`.
De ahí dos predicciones, ambas confirmadas empíricamente:

- **E[v] ≈ 2.0** — medido 1.994–2.002 según la corrida.
- **drift geométrico ≈ 0.75 = 3/4** — medido 0.749–0.751.
- **marginal de `v` ≈ 2^-k** — 0.499, 0.246, 0.124, 0.073, … frente a 0.5, 0.25, 0.125, …

La parte del modelo que **no** es verificable a priori es la **independencia entre pasos**.
Ahí apuntó la investigación.

## 2. La hipótesis

La iteración 1 detectó una **autocorrelación de lag-1 positiva** (≈ +0.0137 en bruto) en
la secuencia de `v` de las trayectorias. Hipótesis: es una **correlación de corto alcance
real**, estructura que la heurística i.i.d. ignora.

## 3. La falsación (iteración 2)

Se montó una prueba falsable con cuatro filtros y una **regla pre-registrada** (ver
`docs/decisions/0001-cierre-veta-lag1.md`):

- **Corrección de sesgo de muestra finita.** Bajo independencia, la ACF muestral no vale 0
  sino `E[ρ_k] ≈ -1/(L-1)`. No se comparó contra 0, sino contra esa línea base. El lag-1
  **corregido** dio **+0.0253** (z ≈ +13.9).
- **Prueba de permutación** (la principal, no asume colas livianas): permutar los `v` dentro
  de cada trayectoria genera el nulo empírico. lag-1 observado significativo, p = 0.0010.
- **Estratificación por longitud `L`** (descartar paradoja de Simpson): el lag-1 corregido se
  mantuvo positivo dentro de los estratos (p.ej. +0.016 en L∈[50,100), +0.031 en [100,200)).
- **Escalamiento por bits:** aquí cayó. El lag-1 corregido **decae monótonamente con la
  magnitud**: +0.052 (20-24 bits) → +0.013 (60-64 bits).

Veredicto provisional: falla la condición 4 → cerrar. **Pero** ese cuarto test tenía un
problema: más bits ⇒ trayectorias más largas, así que **magnitud y longitud `L` estaban
confundidas**. Había que des-confundir antes de ratificar.

## 4. El des-confound (el giro clave)

Este fue el corazón del cierre, y conviene contar también el error de método que apareció
y se corrigió, porque es parte de la honestidad del registro:

- **Confounder de ventana ancha, detectado a mitad de camino.** El primer intento de "fijar
  L" eligió la ventana de L *más ancha* con muestra adecuada — y resultó ser `[11, 261)`,
  que **no fija nada**: cada banda de bits conservaba su L media (de 53 a 187 pasos). Era el
  escalamiento confundido otra vez, disfrazado. Se corrigió el criterio a "ventana **angosta**
  (ancho ≤ 50) que maximiza la muestra", quedando en **[86, 136)** con L media 99–121 entre
  bandas (dispersión 22 frente a ~110 típico): ahí sí, apples-to-apples.

- **Lag-1 a L fija.** Dentro de `[86, 136)`, a **mismo largo**, el lag-1 corregido **sigue
  decayendo con la magnitud**: +0.0511 → +0.0063 (20-24 → 80-84 bits, bandas extendidas).

- **Regresión parcial `lag1 ~ bits + L`** (n = 140,000 trayectorias, des-confound robusto al
  ancho de ventana): coef. de **bits = −0.000624/bit** (IC95 [−0.000659, −0.000589])
  controlando L; coef. de **L = +0.000067/paso** (IC95 > 0). Dos efectos reales y **opuestos**;
  a igual L, la magnitud empuja el lag-1 hacia 0. (Detalle: las bandas altas tienen L media
  *ligeramente mayor* dentro de la ventana, y como el efecto de L es positivo, el confound
  residual empujaría su lag-1 hacia arriba — y aun así quedan más bajas. El efecto real de la
  magnitud es, si acaso, más fuerte de lo que se ve.)

- **Ajuste de asíntota** sobre la serie a L fija: el modelo **"→ 0"** gana por AIC
  (−94.55 vs −92.43); la constante del modelo alternativo es c = −0.0019 (IC95
  [−0.0168, +0.0130], **incluye 0**).

## 5. Veredicto formal

> **La señal lag-1 observada inicialmente no constituye evidencia robusta de una correlación
> estructural persistente en la dinámica impar de Collatz. Tras corregir sesgo de muestra
> finita, controlar por longitud de trayectoria y escalar por magnitud, el efecto decae hacia
> cero. Se cierra como efecto empírico de tamaño finito, no como señal nueva.**

## 6. El patrón multi-lag residual (interpretación, sin experimento nuevo)

El diagnóstico multi-lag mostró que **9/10 lags** sobreviven a la permutación intra-trayectoria
y que el **68% de la energía de la ACF** es dependiente del orden, con el **lag-4 dominando**
al lag-1 (corregido +0.0411 vs +0.0219).

Esto **no** contradice el cierre. La secuencia de `v` es **determinista** dada `n`, así que
correlaciones a varios lags son lo esperable: son la huella de la **mecánica del mapa**, no de
una fuente estocástica nueva. Concretamente, los primeros `k` bits de paridad de una trayectoria
dependen solo de `n mod 2^k` — la **estructura 2-ádica / vector de paridad** ya conocida desde
**Terras (1976)**, *"A stopping time problem on the positive integers"*, Acta Arithmetica 30,
241–252.

**Declaración explícita:** no se persigue una demostración de que el patrón multi-lag es
reducible a "mecánica del mapa + marginal de `v`". Confirmarlo solo ratificaría estructura ya
conocida y **no aportaría descubrimiento**. Queda anotado, no abierto.

## 7. Encuadre: por qué esto cuenta como buen trabajo

Es un **resultado negativo, reproducible y honesto**. El arco —heurística verificada →
hipótesis falsable → regla pre-registrada → cuatro filtros → des-confound (con un error de
método detectado y corregido en el camino)— es exactamente la disciplina que se quería
ejercitar. La maquinaria (ACF agregada con corrección de sesgo, prueba de permutación,
estratificación, regresión parcial de des-confound, ajuste de asíntota por AIC) es
reutilizable en otros sistemas dinámicos.

---

## 8. Congelado para reproducibilidad

**Commit que fija el estado experimental:**
`f57c9abd50ca2eca36d3e2da47acfbcab6989dfc`

**Iteración 2** — `python scripts/run_iteration2.py`
- semilla: `20260601` · B permutaciones: `1000`
- N_eff objetivo banda principal (40-44 bits): `300,000`
- N_eff objetivo por banda de bits: `200,000`
- lags: `1..10` · bandas de bits: `20-24, 30-34, 40-44, 50-54, 60-64`
- min_L: `11`

**Confirmatorio** — `python scripts/run_confirmatorio.py`
- semilla: `20260615` · B permutaciones: `1000`
- trayectorias por banda: `20,000` · umbral N_eff ventana de L: `25,000`
- bandas de bits: `20-24, 30-34, 40-44, 50-54, 60-64, 70-74, 80-84`
- ventana de L fija hallada: `[86, 136)` (ancho máx. permitido 50)
- regresión parcial: `n = 140,000` trayectorias
- supervivencia multi-lag: `n = 5,000` trayectorias (banda 40-44)
- enteros nativos de Python (sin overflow); numba no usado.

**Figuras relevantes** (en `outputs/`, ignorada por git; se regeneran con los scripts y
semillas de arriba — no se regeneró nada para este post-mortem):
- Iteración 1: `trayectorias.png`, `tiempos_parada.png`, `autocorrelacion.png`
- Iteración 2: `it2_acf_bruto.png`, `it2_acf_corregido.png`, `it2_lag1_por_L.png`,
  `it2_lag1_por_bits.png`, `it2_longitudes.png`, `it2_dist_v.png`
- Confirmatorio: `conf_lag1_fijoL.png` (lag-1 a L fija con ajustes de asíntota),
  `conf_multilag.png` (estructura multi-lag que sobrevive a la permutación)

**Reportes completos:** `docs/iteracion2_reporte.md`, `docs/confirmatorio_reporte.md`.
**Decisión formal:** `docs/decisions/0001-cierre-veta-lag1.md`.
