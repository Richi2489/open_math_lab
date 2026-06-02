"""
test_acf.py — pruebas de las herramientas de autocorrelación y significancia
============================================================================

Validan: la ACF por trayectoria, que la corrección de sesgo deja a datos i.i.d.
dentro de banda, que el nulo de permutación se centra en -1/(L-1), y que la
permutación SÍ detecta una autocorrelación positiva inyectada.
"""

import numpy as np
import pytest

from collatz import acf, engine


# ---------------------------------------------------------------------------
# acf_trayectoria
# ---------------------------------------------------------------------------
def test_acf_constante_es_cero():
    a = acf.acf_trayectoria([5, 5, 5, 5], max_lag=2)
    assert a[0] == 0.0  # serie constante: sin variabilidad


def test_acf_alternante_lag1_negativo():
    a = acf.acf_trayectoria([1, -1, 1, -1, 1, -1], max_lag=1)
    assert a[0] == 1.0
    assert a[1] < 0  # señal perfectamente alternante => lag-1 negativo


# ---------------------------------------------------------------------------
# Datos i.i.d. (modelo geométrico P(v=k)=2^-k): la corrección debe dejar ~0
# ---------------------------------------------------------------------------
def _muestra_iid(n_seq, L, rng):
    return [rng.geometric(0.5, size=L).astype(np.float64) for _ in range(n_seq)]


def test_correccion_deja_iid_en_banda():
    rng = np.random.default_rng(123)
    secs = _muestra_iid(400, 40, rng)
    agg = acf.acf_agregado(secs, max_lag=5)
    # El bruto debe ser claramente negativo (sesgo -1/(L-1) ~ -0.026)...
    assert agg["rho_bruto"][1] < 0
    # ...y el corregido debe quedar dentro de unas pocas sigmas de 0.
    assert abs(agg["z"][1]) < 4.0


def test_nulo_permutacion_se_centra_en_baseline():
    rng = np.random.default_rng(7)
    L = 40
    secs = _muestra_iid(300, L, rng)
    res = acf.prueba_permutacion(secs, [1], B=400, rng=np.random.default_rng(8))[1]
    baseline = -1.0 / (L - 1)
    assert abs(res["nulo_media"] - baseline) < 0.02


def test_iid_no_es_significativo_casi_siempre():
    rng = np.random.default_rng(42)
    secs = _muestra_iid(300, 40, rng)
    res = acf.prueba_permutacion(secs, [1], B=400, rng=np.random.default_rng(43))[1]
    # Con datos i.i.d. el observado debe caer dentro del intervalo nulo.
    assert res["p_dos_colas"] > 0.02


# ---------------------------------------------------------------------------
# Autocorrelación positiva inyectada: la permutación DEBE detectarla
# ---------------------------------------------------------------------------
def test_permutacion_detecta_lag1_positivo():
    rng = np.random.default_rng(99)
    secs = []
    for _ in range(200):
        base = rng.geometric(0.5, size=20).astype(np.float64)
        secs.append(np.repeat(base, 2))  # cada valor repetido => fuerte lag-1 positivo
    res = acf.prueba_permutacion(secs, [1], B=400, rng=np.random.default_rng(100))[1]
    assert res["observado"] > res["nulo_media"]
    assert res["significativo"]
    assert res["p_una_cola_sup"] < 0.01


# ---------------------------------------------------------------------------
# Guardas del motor
# ---------------------------------------------------------------------------
def test_overflow_vectorizado_aborta():
    grande = np.array([1 << 62], dtype=np.int64)
    with pytest.raises(OverflowError):
        engine.paso_acelerado_vector(grande)


def test_max_iter_dispara():
    with pytest.raises(RuntimeError):
        engine.secuencia_v(27, max_iter=3)
