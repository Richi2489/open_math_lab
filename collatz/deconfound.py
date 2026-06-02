"""
deconfound.py — des-confundir magnitud (bits) de largo (L)
==========================================================

El experimento confirmatorio de la iteración 2 mostró que el lag-1 corregido decae
con los bits del arranque. Pero bits y largo L están correlacionados (más bits =>
trayectorias más largas), así que ese decaimiento podía deberse a L, no a la magnitud.

Aquí se aísla un efecto del otro de dos formas independientes:

1. VENTANA DE L FIJA (`buscar_ventana_L`, `lag1_por_banda_en_ventana`): se elige la
   ventana de L más ancha con muestra adecuada en TODAS las bandas de bits y, dentro
   de ella, se recalcula el lag-1 por banda. Misma L => comparación apples-to-apples.

2. REGRESIÓN PARCIAL (`regresion_parcial`): regresión lineal ponderada del lag-1
   corregido por trayectoria sobre (bits, L) usando TODOS los datos. El coeficiente
   de `bits` es el efecto de la magnitud CONTROLANDO por L. Robusto al ancho de la
   ventana.

Más:
  - `ajustar_asintota`: ajusta lag-1(bits) a un decaimiento a 0 vs a una constante c>0.
  - `supervivencia_multilag`: fracción del patrón ACF (lags 1..10) que sobrevive a la
    permutación intra-trayectoria (= estructura dependiente del orden = mecánica).
"""

from __future__ import annotations

import numpy as np

from . import acf


def _longitudes(seqs) -> np.ndarray:
    return np.array([len(s) for s in seqs], dtype=np.int64)


# ---------------------------------------------------------------------------
# 1. Ventana de L común a todas las bandas
# ---------------------------------------------------------------------------
def buscar_ventana_L(seqs_por_banda: dict, umbral_neff: int, max_ancho: int = 50,
                     lo_min: int = 11, hi_max: int = 300, paso: int = 5):
    """Encuentra una ventana [lo, hi) de L que CONTROLE el largo (ancho acotado) y
    maximice la muestra mínima entre bandas.

    Criterio: entre todas las ventanas con ancho <= `max_ancho` cuya N_eff(lag1)
    mínima sobre las bandas sea >= `umbral_neff`, devuelve la que MAXIMIZA esa N_eff
    mínima. Limitar el ancho es clave: una ventana ancha (p.ej. [11, 261)) no fija L
    —cada banda conserva su L media—, así que no des-confunde. Una ventana angosta en
    la zona de solape sí deja L aproximadamente constante entre bandas.

    Devuelve (lo, hi) o None si ninguna ventana acotada cumple el umbral en todas.
    """
    Ls = {b: _longitudes(s) for b, s in seqs_por_banda.items()}
    bordes = list(range(lo_min, hi_max + paso, paso))
    mejor = None  # (min_neff, lo, hi)
    for i in range(len(bordes)):
        for j in range(i + 1, len(bordes)):
            lo, hi = bordes[i], bordes[j]
            if hi - lo > max_ancho:
                break
            min_neff = np.inf
            for L in Ls.values():
                m = (L >= lo) & (L < hi)
                neff = int((L[m] - 1).sum()) if m.any() else 0
                if neff < min_neff:
                    min_neff = neff
                if min_neff < umbral_neff:
                    break
            if min_neff >= umbral_neff and (mejor is None or min_neff > mejor[0]):
                mejor = (min_neff, lo, hi)
    if mejor is None:
        return None
    return mejor[1], mejor[2]


def subconjunto_en_ventana(seqs, lo: int, hi: int):
    return [s for s in seqs if lo <= len(s) < hi]


def lag1_por_banda_en_ventana(seqs_por_banda: dict, lo: int, hi: int,
                              B: int, rng) -> dict:
    """Lag-1 corregido + permutación por banda, restringido a L en [lo, hi).

    Devuelve {banda: dict(n, L_medio, n_eff, lag1_corregido, banda_signif, perm)}.
    Se reporta L_medio por banda como diagnóstico de confound RESIDUAL: si difiere
    mucho entre bandas, la ventana es demasiado ancha y conviene leer la regresión.
    """
    res = {}
    for b, seqs in seqs_por_banda.items():
        sub = subconjunto_en_ventana(seqs, lo, hi)
        agg = acf.acf_agregado(sub, max_lag=1)
        perm = acf.prueba_permutacion(sub, [1], B=B, rng=rng)[1]
        Ls = _longitudes(sub)
        res[b] = {
            "n": len(sub),
            "L_medio": float(Ls.mean()) if len(sub) else float("nan"),
            "n_eff": int(agg["n_eff"][1]),
            "lag1_corregido": float(agg["rho_corregido"][1]),
            "banda_signif": float(agg["banda"][1]),
            "perm": perm,
        }
    return res


# ---------------------------------------------------------------------------
# 2. Regresión parcial: efecto de bits controlando por L
# ---------------------------------------------------------------------------
def regresion_parcial(seqs_por_banda: dict, centros_bits: dict) -> dict:
    """Regresión lineal ponderada del lag-1 corregido por trayectoria sobre bits y L.

    Modelo:  corr_i = b0 + b_bits * bits_i + b_L * L_i + ruido,   peso_i = L_i - 1.
    corr_i = rho1_i + 1/(L_i - 1)  (lag-1 con el sesgo de muestra finita descontado).

    Devuelve dict con beta, errores estándar (se), IC 95% y n. El coeficiente de
    interés es b_bits: efecto de la magnitud A IGUAL L.
    """
    ys, b_bits, b_L, w = [], [], [], []
    for b, seqs in seqs_por_banda.items():
        c = centros_bits[b]
        for s in seqs:
            L = len(s)
            if L < 3:
                continue
            rho1 = acf.acf_trayectoria(s, 1)[1]
            ys.append(rho1 + 1.0 / (L - 1))
            b_bits.append(c)
            b_L.append(L)
            w.append(L - 1)
    ys = np.asarray(ys, dtype=np.float64)
    X = np.column_stack([np.ones(len(ys)), np.asarray(b_bits, float), np.asarray(b_L, float)])
    w = np.asarray(w, dtype=np.float64)
    sw = np.sqrt(w)
    Xw, yw = X * sw[:, None], ys * sw
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    resid = yw - Xw @ beta
    dof = len(ys) - X.shape[1]
    sigma2 = float(resid @ resid) / dof
    cov = sigma2 * np.linalg.inv(Xw.T @ Xw)
    se = np.sqrt(np.diag(cov))
    return {
        "nombres": ["intercepto", "bits", "L"],
        "beta": beta,
        "se": se,
        "ic95": np.column_stack([beta - 1.96 * se, beta + 1.96 * se]),
        "n": len(ys),
    }


# ---------------------------------------------------------------------------
# 3. Ajuste de asíntota: decaimiento a 0 vs a constante c>0
# ---------------------------------------------------------------------------
def _modelo_cero(x, a, k):
    return a * np.exp(-k * x)


def _modelo_const(x, c, a, k):
    return c + a * np.exp(-k * x)


def ajustar_asintota(bits, y, sigma=None) -> dict:
    """Ajusta lag-1(bits) a dos modelos y los compara por AIC:
      (a) a·exp(-k·bits)      -> asíntota 0
      (b) c + a·exp(-k·bits)  -> asíntota c

    Devuelve dict con params, RSS, AIC de cada modelo, el ganador y la asíntota c
    (modelo b) con su error estándar. Requiere scipy.
    """
    from scipy.optimize import curve_fit

    x = np.asarray(bits, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    x0 = x - x.min()  # estabiliza el ajuste
    n = len(x)
    abs_sigma = sigma is not None

    def _aic(rss, p):
        return n * np.log(rss / n) + 2 * p

    out = {"x_min": float(x.min())}
    try:
        pa, _ = curve_fit(_modelo_cero, x0, y, p0=[y[0], 0.05],
                          sigma=sigma, absolute_sigma=abs_sigma, maxfev=20000)
        rss_a = float(np.sum((y - _modelo_cero(x0, *pa)) ** 2))
        out["modelo_cero"] = {"params": pa, "rss": rss_a, "aic": _aic(rss_a, 2),
                              "asintota": 0.0}
    except Exception as e:  # pragma: no cover
        out["modelo_cero"] = {"error": str(e), "aic": np.inf}

    try:
        pb, cb = curve_fit(_modelo_const, x0, y,
                           p0=[y[-1], y[0] - y[-1], 0.05],
                           sigma=sigma, absolute_sigma=abs_sigma, maxfev=20000)
        rss_b = float(np.sum((y - _modelo_const(x0, *pb)) ** 2))
        se_c = float(np.sqrt(cb[0, 0])) if np.all(np.isfinite(cb)) else float("nan")
        out["modelo_const"] = {"params": pb, "rss": rss_b, "aic": _aic(rss_b, 3),
                               "asintota": float(pb[0]), "se_asintota": se_c}
    except Exception as e:  # pragma: no cover
        out["modelo_const"] = {"error": str(e), "aic": np.inf}

    out["ganador"] = "cero" if out["modelo_cero"]["aic"] <= out["modelo_const"]["aic"] else "const"
    return out


# ---------------------------------------------------------------------------
# 4. Supervivencia multi-lag a la permutación
# ---------------------------------------------------------------------------
def supervivencia_multilag(secuencias, lags, B: int, rng) -> dict:
    """Cuantifica cuánto del patrón ACF (lags 1..K) sobrevive a la permutación
    intra-trayectoria (= estructura dependiente del ORDEN = mecánica del mapa).

    Por lag: observado, nulo medio, "corregido" = observado - nulo medio, significancia.
    Agregado: energía observada sum_k obs_k^2 vs distribución nula de energía;
    fracción de energía que sobrevive y p-valor empírico.
    """
    perm = acf.prueba_permutacion(secuencias, lags, B=B, rng=rng, devolver_nulo=True)
    nulos = perm["_nulo"]
    por_lag = {}
    e_obs = 0.0
    e_null = np.zeros(B)
    n_sig = 0
    for k in lags:
        d = perm[k]
        por_lag[k] = {
            "observado": d["observado"],
            "nulo_media": d["nulo_media"],
            "corregido": d["observado"] - d["nulo_media"],
            "significativo": d["significativo"],
            "p_dos_colas": d["p_dos_colas"],
        }
        n_sig += int(d["significativo"])
        e_obs += d["observado"] ** 2
        e_null += nulos[k] ** 2
    frac = float((e_obs - e_null.mean()) / e_obs) if e_obs > 0 else 0.0
    p_energia = (1 + int(np.sum(e_null >= e_obs))) / (B + 1)
    return {
        "por_lag": por_lag,
        "energia_obs": float(e_obs),
        "energia_nulo_media": float(e_null.mean()),
        "fraccion_energia_sobrevive": frac,
        "p_energia": p_energia,
        "lags_significativos": n_sig,
        "n_lags": len(lags),
    }
