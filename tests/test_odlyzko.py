"""
test_odlyzko.py — pruebas del loader de Odlyzko y la tubería de alta altura
===========================================================================

No tocan la ruta de mpmath (`riemann.zeros`); usan archivos temporales y un fixture
pequeño y realista (primeros 500 ceros reales) en tests/fixtures/.
"""

from pathlib import Path

import numpy as np
import pytest

from riemann import gue, metrics, odlyzko, spacing

_FIXTURE = Path(__file__).parent / "fixtures" / "odlyzko_small.txt"


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------
def test_odlyzko_loader_txt_one_column(tmp_path):
    p = tmp_path / "z.txt"
    p.write_text("   14.134725142\n   21.022039639\n   25.010857580\n", encoding="utf-8")
    g = odlyzko.cargar_odlyzko(p)
    assert g.dtype == np.float64 and g.ndim == 1
    assert abs(g[0] - 14.134725142) < 1e-9
    assert len(g) == 3


def test_odlyzko_loader_csv_named_column(tmp_path):
    p = tmp_path / "z.csv"
    p.write_text("n,gamma\n1,14.134725\n2,21.022040\n3,25.010858\n", encoding="utf-8")
    g = odlyzko.cargar_odlyzko(p)
    assert abs(g[0] - 14.134725) < 1e-6
    assert len(g) == 3


def test_odlyzko_rejects_unsorted():
    with pytest.raises(ValueError):
        odlyzko.validar_zeros(np.array([14.0, 21.0, 20.0, 25.0]))


def test_odlyzko_rejects_duplicates():
    with pytest.raises(ValueError):
        odlyzko.validar_zeros(np.array([14.0, 21.0, 21.0, 25.0]))


def test_odlyzko_rejects_nan():
    with pytest.raises(ValueError):
        odlyzko.validar_zeros(np.array([14.0, np.nan, 25.0]))


def test_odlyzko_summary_height_range():
    g = np.array([14.0, 21.0, 25.0, 30.0])
    r = odlyzko.resumir_dataset(g)
    assert r["n"] == 4
    assert r["T_min"] == 14.0 and r["T_max"] == 30.0


# ---------------------------------------------------------------------------
# Unfolding y ventanas
# ---------------------------------------------------------------------------
def test_local_density_unfold_positive_mean_near_one():
    g = odlyzko.cargar_odlyzko(_FIXTURE)
    s = spacing.unfold_densidad_local(g)
    assert np.all(s > 0)
    assert abs(s.mean() - 1.0) < 0.05


def test_window_split_nonoverlapping():
    g = np.arange(100.0)
    ventanas = odlyzko.dividir_en_ventanas(g, 4)
    assert len(ventanas) == 4
    # reconstrucción exacta y sin solape
    assert np.array_equal(np.concatenate(ventanas), g)
    assert sum(len(v) for v in ventanas) == len(g)


def test_odlyzko_pipeline_small_fixture():
    g = odlyzko.cargar_odlyzko(_FIXTURE)
    odlyzko.validar_zeros(g)
    s = spacing.unfold_densidad_local(g)
    rng = np.random.default_rng(0)
    s_gue = gue.espaciados_gue(120, rng, n_matrices=10)
    d = metrics.distancias(s, s_gue)
    assert abs(s.mean() - 1.0) < 0.05
    assert 0.0 <= d["ks"] <= 1.0
    assert np.isfinite(d["wasserstein"])
