import numpy as np
import pytest
from synth import stefan_columns, truncate

from exocam_accelerate.advise import (TAPERED_ADVICE_SCHEMA_VERSION, AdvisorConfig,
                                      advise, taper_advice)
from exocam_accelerate.aqua_ice import AquaIcePlugin
from exocam_accelerate.check import Verdict, check_jump
from exocam_accelerate.taper import (TaperConfig, cell_factors, effective_factor,
                                     imbalance_ratio, solve_peak, stefan_mask)

NJ, NI = 4, 6
DT = 10.0


def grid(frac_static=0.0, partial=0, h_static=10.0):
    """Stefan-growing ice (S = 20 m^2/yr) with the first columns static
    (near local equilibrium) and ``partial`` cells of partial ice cover."""
    h_old = np.full((NJ, NI), 50.0)
    h_now = np.sqrt(h_old**2 + 2 * 20.0 * DT)
    n_static = int(round(frac_static * NI))
    h_old[:, :n_static] = h_now[:, :n_static] = h_static
    aice = np.ones((NJ, NI))
    aice.flat[NJ * NI - partial:] = 0.5 if partial else 1.0
    area = np.ones((NJ, NI))
    return h_old, h_now, aice, area


def test_stefan_growing_ice_gets_full_weight_static_gets_none():
    h_old, h_now, aice, area = grid(1 / 3)
    m = stefan_mask(h_old, h_now, DT, aice, area)
    assert np.all(m.weight[:, 2:] == 1.0)
    assert np.all(m.weight[:, :2] == 0.0)
    s = m.summary()
    assert s["area_full_weight"] == pytest.approx(2 / 3)
    assert s["area_unscaled"] == pytest.approx(1 / 3)


def test_partial_ice_cover_is_never_scaled():
    h_old, h_now, aice, area = grid(partial=3)
    m = stefan_mask(h_old, h_now, DT, aice, area)
    assert np.all(m.weight.flat[-3:] == 0.0)
    assert np.all(m.weight.flat[:-3] == 1.0)


def test_ramp_is_linear_between_lo_and_hi():
    h_old, h_now, aice, area = grid()
    # halve the Stefan product in one cell, quarter it in another
    for (j, i), s in (((0, 0), 0.5), ((0, 1), 0.1)):
        h_now[j, i] = np.sqrt(h_old[j, i] ** 2 + 2 * 20.0 * DT * s)
    m = stefan_mask(h_old, h_now, DT, aice, area, TaperConfig(0.0, 1.0))
    assert m.w[0, 0] == pytest.approx(0.5, rel=1e-2)
    assert m.weight[0, 0] == pytest.approx(m.w[0, 0])
    m2 = stefan_mask(h_old, h_now, DT, aice, area, TaperConfig(0.2, 0.4))
    assert m2.weight[0, 1] == 0.0 and m2.weight[0, 0] == 1.0


def test_refuses_when_thick_ice_not_growing():
    h = np.full((NJ, NI), 30.0)
    with pytest.raises(ValueError, match="not growing"):
        stefan_mask(h, h, DT, np.ones_like(h), np.ones_like(h))


@pytest.mark.parametrize("kw", [dict(ramp_lo=0.5, ramp_hi=0.4),
                                dict(ramp_lo=-0.1), dict(full_cover=0.0)])
def test_invalid_config(kw):
    with pytest.raises(ValueError):
        TaperConfig(**kw)


def test_uniform_weight_reduces_to_global_law():
    m = stefan_mask(*grid()[:2], DT, *grid()[2:])
    for F in (1.2, 1.5):
        assert imbalance_ratio(F, m) == pytest.approx(1 / F)
        assert effective_factor(F, m) == pytest.approx(F)
        assert solve_peak(1 / F, m, 2.0) == pytest.approx(F, rel=1e-8)


def test_partial_taper_needs_a_larger_peak():
    h_old, h_now, aice, area = grid(1 / 3)
    m = stefan_mask(h_old, h_now, DT, aice, area)
    # static cells carry no freezing: the growing ones alone set N, so R = 1/F
    assert imbalance_ratio(1.5, m) == pytest.approx(1 / 1.5)
    assert effective_factor(1.5, m) < 1.5
    np.testing.assert_allclose(cell_factors(1.5, m.weight)[:, :2], 1.0)


def test_unreachable_target_returns_none():
    h_old, h_now, aice, area = grid(1 / 3)
    h_now[:, :2] = h_old[:, :2] + 5.0           # static-ish but growing a bit
    m = stefan_mask(h_old, h_now, DT, aice, area, TaperConfig(0.9, 1.0))
    assert solve_peak(1e-3, m, 2.0) is None


# ---- advice ----

def uniform_advice():
    adv = advise(truncate(stefan_columns(), 60), "synth", AdvisorConfig(max_ice_factor=1.5))
    assert adv.jump, adv.reasons
    return adv.to_dict()


def test_taper_advice_with_uniform_mask_matches_uniform_advice():
    u = uniform_advice()
    m = stefan_mask(*grid()[:2], DT, *grid()[2:])
    t = taper_advice(u, m)
    assert t["schema_version"] == TAPERED_ADVICE_SCHEMA_VERSION
    assert t["ice_factor"] == pytest.approx(u["ice_factor"])
    assert t["N_after"] == pytest.approx(u["N_after"])
    assert t["effective_factor"] == pytest.approx(u["ice_factor"])
    assert t["law_offset_after"] == pytest.approx(0.0, abs=1e-9)
    assert t["years_skipped"] == pytest.approx(u["years_skipped"])
    assert t["expected_after_jump"]["TS"] == pytest.approx(u["expected_after_jump"]["TS"])
    assert t["taper"]["uniform"]["N_after"] == u["N_after"]


def test_taper_advice_partial_mask():
    u = uniform_advice()
    m = stefan_mask(*grid(1 / 3)[:2], DT, *grid(1 / 3)[2:])
    t = taper_advice(u, m)
    assert t["effective_factor"] < t["ice_factor"]
    # thinner global-mean ice than uniform, same freezing-carried N: the run
    # sits above (less negative than) the global law at that thickness
    assert t["law_offset_after"] > 0
    assert t["N_after"] == pytest.approx(u["N_after"])


def test_taper_advice_refusals():
    m = stefan_mask(*grid()[:2], DT, *grid()[2:])
    with pytest.raises(ValueError, match="do not jump"):
        taper_advice({"ice_factor": None}, m)
    t = taper_advice(uniform_advice(), m)
    with pytest.raises(ValueError, match="already tapered"):
        taper_advice(t, m)


# ---- plugin and check ----

def test_plugin_applies_factor_map_per_cell():
    rng = np.random.default_rng(0)
    f = {"vicen": rng.uniform(0, 30, (5, NJ, NI)),
         "eicen": rng.uniform(-1e9, -1e8, (20, NJ, NI))}
    fmap = np.ones((NJ, NI))
    fmap[:, 3:] = 1.5
    out = AquaIcePlugin().apply_delta(f, (fmap, 1.0))
    np.testing.assert_allclose(out["vicen"], f["vicen"] * fmap)
    np.testing.assert_allclose(out["eicen"] / out["vicen"].repeat(4, axis=0),
                               f["eicen"] / f["vicen"].repeat(4, axis=0))
    with pytest.raises(ValueError, match="hard bound"):
        AquaIcePlugin().apply_delta(f, (fmap * 2, 1.0))


def test_check_uses_effective_factor_and_offset():
    u = uniform_advice()
    u.update(effective_factor=1.3, law_offset_after=1.0)
    log = {"jump_model_year": 61, "ice_factor": 1.5, "advice": u}
    cols = stefan_columns(years=70, jump_year=61, factor=1.3, N_offset_after=1.0)
    res = check_jump(cols, log)
    assert res.verdict is Verdict.PASS, res.reasons
    assert abs(res.metrics["land_error"]) < 0.02
    # the same run judged as a uniform x1.5 jump would not pass
    plain = dict(log, advice=dict(u, effective_factor=None, law_offset_after=None))
    assert check_jump(cols, plain).verdict is Verdict.FAIL
