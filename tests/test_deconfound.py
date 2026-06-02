"""
test_deconfound.py — pruebas de la des-confusión bits vs L
==========================================================
"""

import numpy as np
import pytest

from collatz import deconfound


def _seqs_de_longitudes(longitudes):
    """Secuencias dummy con las longitudes dadas (el contenido no importa aquí)."""
    return [np.ones(int(L)) for L in longitudes]


# ---------------------------------------------------------------------------
# Búsqueda de ventana de L
# ---------------------------------------------------------------------------
def test_buscar_ventana_encuentra_solape():
    rng = np.random.default_rng(0)
    bandas = {
        "baja": _seqs_de_longitudes(rng.integers(40, 120, 4000)),
        "media": _seqs_de_longitudes(rng.integers(80, 160, 4000)),
        "alta": _seqs_de_longitudes(rng.integers(100, 200, 4000)),
    }
    v = deconfound.buscar_ventana_L(bandas, umbral_neff=5000, paso=5)
    assert v is not None
    lo, hi = v
    assert 0 < lo < hi


def test_buscar_ventana_umbral_inalcanzable_devuelve_none():
    # La banda corta aporta como máximo 500*(20-1)=9500 pares: ninguna ventana
    # puede darle N_eff>=20000, así que no hay ventana común.
    bandas = {
        "baja": _seqs_de_longitudes([20] * 500),
        "alta": _seqs_de_longitudes([200] * 5000),
    }
    assert deconfound.buscar_ventana_L(bandas, umbral_neff=20000, paso=5) is None


# ---------------------------------------------------------------------------
# Regresión parcial
# ---------------------------------------------------------------------------
def _seq_iid(L, rng):
    return rng.geometric(0.5, size=L).astype(np.float64)


def _seq_bloques(L_pares, rng):
    base = rng.geometric(0.5, size=L_pares).astype(np.float64)
    return np.repeat(base, 2)  # lag-1 fuertemente positivo


def test_regresion_detecta_efecto_bits():
    rng = np.random.default_rng(1)
    seqs = {}
    centros = {}
    # banda de bits baja: i.i.d. (lag-1 ~ 0); banda alta: bloques (lag-1 alto).
    seqs["10"] = [_seq_iid(int(rng.integers(30, 50)), rng) for _ in range(400)]
    seqs["30"] = [_seq_bloques(int(rng.integers(15, 25)), rng) for _ in range(400)]
    centros = {"10": 10.0, "30": 30.0}
    reg = deconfound.regresion_parcial(seqs, centros)
    # coef. de bits debe ser positivo y significativo (IC95 enteramente > 0).
    assert reg["ic95"][1][0] > 0


# ---------------------------------------------------------------------------
# Ajuste de asíntota
# ---------------------------------------------------------------------------
def test_asintota_constante_positiva():
    x = np.array([20, 30, 40, 50, 60, 70, 80], dtype=float)
    y = 0.05 + 0.10 * np.exp(-0.10 * (x - 20))
    aj = deconfound.ajustar_asintota(x, y)
    assert abs(aj["modelo_const"]["asintota"] - 0.05) < 0.02


def test_asintota_hacia_cero():
    x = np.array([20, 30, 40, 50, 60, 70, 80], dtype=float)
    y = 0.20 * np.exp(-0.08 * (x - 20))
    aj = deconfound.ajustar_asintota(x, y)
    # la constante ajustada debe quedar cerca de 0.
    assert abs(aj["modelo_const"]["asintota"]) < 0.03


# ---------------------------------------------------------------------------
# Supervivencia multi-lag
# ---------------------------------------------------------------------------
def test_supervivencia_alta_en_datos_ordenados():
    rng = np.random.default_rng(2)
    seqs = [_seq_bloques(15, rng) for _ in range(150)]
    res = deconfound.supervivencia_multilag(
        seqs, [1, 2, 3], B=200, rng=np.random.default_rng(3)
    )
    assert res["fraccion_energia_sobrevive"] > 0.5
    assert res["p_energia"] < 0.05


def test_supervivencia_baja_en_iid():
    rng = np.random.default_rng(4)
    seqs = [_seq_iid(40, rng) for _ in range(200)]
    res = deconfound.supervivencia_multilag(
        seqs, [1, 2, 3], B=200, rng=np.random.default_rng(5)
    )
    # En i.i.d. la energía observada no debe exceder mucho a la nula.
    assert res["p_energia"] > 0.02
