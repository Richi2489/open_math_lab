"""
acf.py — autocorrelación de los v's y su significancia
======================================================

Herramientas para la iteración 2: decidir, con rigor falsable, si la
autocorrelación de lag-1 en la secuencia de v's del mapa acelerado es señal real
o artefacto de muestra finita.

Definiciones y convenciones
---------------------------
ACF muestral por trayectoria (lag k):

    rho_k = sum_{t} (v_t - vbar)(v_{t+k} - vbar) / sum_t (v_t - vbar)^2

con vbar la media de ESA trayectoria.

Agregado entre trayectorias: media de los rho_k por trayectoria PONDERADA por el
número de pares disponibles en ese lag, w_i = max(L_i - k, 0). Esto hace que
N_eff(k) = sum_i w_i sea el número total de pares y la banda analítica
±1.96/sqrt(N_eff) tenga el significado usual.

Sesgo de muestra finita: bajo independencia, E[rho_k] ≈ -1/(L-1) por trayectoria
(no 0). NO se compara contra cero: se compara contra esa línea base. El ACF
"corregido" resta -1/(L_i-1) a cada trayectoria antes de agregar.

CAVEAT: la banda analítica supone colas bien portadas. Los v's tienen cola pesada
(geométrica), así que la prueba PRINCIPAL de significancia es de PERMUTACIÓN
(`prueba_permutacion`), que no asume nada sobre las colas.
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# ACF de una trayectoria
# ---------------------------------------------------------------------------
def acf_trayectoria(v, max_lag: int) -> np.ndarray:
    """ACF muestral (lags 0..max_lag) de una secuencia 1D. acf[0] = 1 (o 0 si es
    constante). Los lags sin soporte (k >= L) quedan en 0."""
    x = np.asarray(v, dtype=np.float64)
    L = len(x)
    x = x - x.mean()
    denom = float(np.dot(x, x))
    acf = np.zeros(max_lag + 1)
    if denom == 0.0:
        return acf  # serie constante: sin variabilidad
    acf[0] = 1.0
    for k in range(1, min(max_lag, L - 1) + 1):
        acf[k] = float(np.dot(x[:-k], x[k:])) / denom
    return acf


# ---------------------------------------------------------------------------
# ACF agregado con corrección de sesgo y banda analítica
# ---------------------------------------------------------------------------
def acf_agregado(secuencias, max_lag: int) -> dict:
    """Agrega la ACF sobre muchas trayectorias.

    Devuelve un dict con arrays indexados por lag (0..max_lag):
        lags, rho_bruto, rho_corregido, baseline_medio, n_eff, banda, z
    donde:
        rho_bruto[k]     : media ponderada de rho_k por trayectoria.
        baseline_medio[k]: media ponderada de -1/(L_i-1).
        rho_corregido[k] : media ponderada de (rho_k + 1/(L_i-1)).
        n_eff[k]         : sum_i max(L_i-k, 0)  (número de pares).
        banda[k]         : 1.96/sqrt(n_eff[k])  (analítica, ~gaussiana).
        z[k]             : rho_corregido[k] * sqrt(n_eff[k]).
    """
    num_bruto = np.zeros(max_lag + 1)
    num_corr = np.zeros(max_lag + 1)
    base_acc = np.zeros(max_lag + 1)
    n_eff = np.zeros(max_lag + 1)

    for v in secuencias:
        x = np.asarray(v, dtype=np.float64)
        L = len(x)
        if L < 2:
            continue
        x = x - x.mean()
        denom = float(np.dot(x, x))
        if denom == 0.0:
            continue
        baseline = -1.0 / (L - 1)
        for k in range(1, max_lag + 1):
            w = L - k
            if w <= 0:
                continue
            rho_k = float(np.dot(x[:-k], x[k:])) / denom
            num_bruto[k] += w * rho_k
            num_corr[k] += w * (rho_k - baseline)
            base_acc[k] += w * baseline
            n_eff[k] += w

    with np.errstate(invalid="ignore", divide="ignore"):
        rho_bruto = np.where(n_eff > 0, num_bruto / n_eff, 0.0)
        rho_corr = np.where(n_eff > 0, num_corr / n_eff, 0.0)
        base_medio = np.where(n_eff > 0, base_acc / n_eff, 0.0)
        banda = np.where(n_eff > 0, 1.96 / np.sqrt(n_eff), np.nan)
        z = np.where(n_eff > 0, rho_corr * np.sqrt(n_eff), 0.0)

    rho_bruto[0] = 1.0
    rho_corr[0] = 0.0
    return {
        "lags": np.arange(max_lag + 1),
        "rho_bruto": rho_bruto,
        "rho_corregido": rho_corr,
        "baseline_medio": base_medio,
        "n_eff": n_eff,
        "banda": banda,
        "z": z,
    }


# ---------------------------------------------------------------------------
# Prueba de permutación (PRINCIPAL)
# ---------------------------------------------------------------------------
def prueba_permutacion(secuencias, lags, B: int = 1000, rng=None,
                       devolver_nulo: bool = False) -> dict:
    """Distribución nula empírica del ACF agregado por permutación.

    Para cada trayectoria se permuta aleatoriamente su secuencia de v's (lo que
    destruye cualquier orden temporal pero conserva la distribución marginal de v
    y la longitud L, así que el nulo incorpora AUTOMÁTICAMENTE el sesgo -1/(L-1) y
    la cola pesada). Se repite B veces y se agrega igual que el observado.

    Devuelve {k: dict(...)} con, por cada lag k:
        observado     : rho_k agregado (bruto) en los datos reales.
        nulo_media    : media de la distribución nula del rho_k agregado.
        p2_5, p97_5   : percentiles 2.5 y 97.5 del nulo.
        p_dos_colas   : p-valor empírico de dos colas (centrado en nulo_media).
        p_una_cola_sup: p-valor empírico de una cola (H1: rho_k > nulo).
        significativo : observado fuera de [p2_5, p97_5].

    Si `devolver_nulo=True`, añade además la clave especial '_nulo': {k: array(B)}
    con la distribución nula completa del rho_k agregado (para análisis de energía).
    """
    if rng is None:
        rng = np.random.default_rng()
    lags = [int(k) for k in lags]
    num_obs = {k: 0.0 for k in lags}
    peso = {k: 0.0 for k in lags}
    null = {k: np.zeros(B) for k in lags}

    for v in secuencias:
        x = np.asarray(v, dtype=np.float64)
        L = len(x)
        xc = x - x.mean()
        denom = float(np.dot(xc, xc))
        if denom == 0.0:
            continue

        # --- observado ---
        for k in lags:
            w = L - k
            if w <= 0:
                continue
            num_obs[k] += w * (float(np.dot(xc[:-k], xc[k:])) / denom)
            peso[k] += w

        # --- B permutaciones de esta trayectoria (vectorizado sobre B) ---
        M = np.broadcast_to(x, (B, L)).copy()
        M = rng.permuted(M, axis=1)
        Mc = M - M.mean(axis=1, keepdims=True)
        den = np.einsum("ij,ij->i", Mc, Mc)  # (B,)
        seguro = den > 0
        for k in lags:
            w = L - k
            if w <= 0:
                continue
            numk = np.einsum("ij,ij->i", Mc[:, :-k], Mc[:, k:])
            rho_perm = np.zeros(B)
            rho_perm[seguro] = numk[seguro] / den[seguro]
            null[k] += w * rho_perm

    resultados = {}
    nulos_completos = {}
    for k in lags:
        if peso[k] == 0:
            continue
        obs = num_obs[k] / peso[k]
        nulo = null[k] / peso[k]
        nulos_completos[k] = nulo
        media = float(nulo.mean())
        p_lo = float(np.percentile(nulo, 2.5))
        p_hi = float(np.percentile(nulo, 97.5))
        dif = abs(obs - media)
        p_dos = (1 + int(np.sum(np.abs(nulo - media) >= dif))) / (B + 1)
        p_sup = (1 + int(np.sum(nulo >= obs))) / (B + 1)
        resultados[k] = {
            "observado": obs,
            "nulo_media": media,
            "nulo_std": float(nulo.std(ddof=1)),
            "p2_5": p_lo,
            "p97_5": p_hi,
            "p_dos_colas": p_dos,
            "p_una_cola_sup": p_sup,
            "significativo": (obs < p_lo) or (obs > p_hi),
        }
    if devolver_nulo:
        resultados["_nulo"] = nulos_completos
    return resultados
