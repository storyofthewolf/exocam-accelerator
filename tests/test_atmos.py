import json

import numpy as np
import pytest

from exocam_accelerate.atmos import (AtmConstants, CamAtmPlugin, ProfileConfig,
                                     STATE_FIELDS, esat, interfaces, jump_state,
                                     measured_profile, p_mid, potential_temperature_field,
                                     qsat, temperature, total_energy)

C = AtmConstants(cpair=1024.26, rair=288.70, zvir=0.5986, gravit=9.8, ptop=3.26)
NL, NJ, NI = 8, 3, 4


def state(seed=0, qscale=1.0):
    rng = np.random.default_rng(seed)
    ps = 4.3e5 + 1e3 * rng.standard_normal((NJ, NI))
    frac = np.linspace(0.02, 1.0, NL + 1)                  # interface fractions
    pe = C.ptop + (ps - C.ptop)[None] * (frac[:, None, None] ** 2)
    pe[0] = C.ptop
    D = np.diff(pe, axis=0)
    pm = p_mid(pe)
    T = 200.0 + 170.0 * (pm / pm[-1]) ** 0.3 + rng.normal(0, 0.5, pm.shape)
    Q = np.clip(0.5 * qsat(T, pm, C.epsilon) * qscale, 1e-7, 0.4)
    st = {"DELP": D, "Q": Q, "PT": potential_temperature_field(T, Q, D, C),
          "PS": pe[-1], "U": rng.normal(0, 10, D.shape), "V": rng.normal(0, 5, D.shape),
          "PHIS": np.zeros((NJ, NI)), "CLDLIQ": np.full(D.shape, 1e-5),
          "CLDICE": np.full(D.shape, 1e-6), "TCWAT": T + 0.1, "QCWAT": Q.copy(),
          "LCWAT": np.full(D.shape, 1.1e-5), "T_TTEND": T.copy()}
    st["TEOUT"] = total_energy(T, Q, D, st["U"], st["V"], st["PHIS"], st["CLDLIQ"], C)
    return st, T


def test_temperature_round_trip():
    st, T = state()
    np.testing.assert_allclose(temperature(st["PT"], st["Q"], st["DELP"], C), T, rtol=1e-12)


def test_esat_sane():
    assert esat(373.124) == pytest.approx(101325, rel=2e-3)
    assert esat(273.16) == pytest.approx(611.66, rel=5e-3)


def test_jump_state_physics():
    st, T = state()
    dT = np.linspace(0, 2, NL)                             # more warming low down
    out, rep = jump_state(st, dT, C)
    Tn = temperature(out["PT"], out["Q"], out["DELP"], C)
    np.testing.assert_allclose(Tn - T, np.broadcast_to(dT[:, None, None], T.shape), atol=1e-9)
    # dry mass per layer kept; PS consistent with DELP
    np.testing.assert_allclose(out["DELP"] * (1 - out["Q"]), st["DELP"] * (1 - st["Q"]),
                               rtol=1e-12)
    np.testing.assert_allclose(out["PS"], interfaces(out["DELP"], C.ptop)[-1])
    # fixed RH
    pm0, pm1 = p_mid(interfaces(st["DELP"], C.ptop)), p_mid(interfaces(out["DELP"], C.ptop))
    rh0 = st["Q"] / qsat(T, pm0, C.epsilon)
    rh1 = out["Q"] / qsat(Tn, pm1, C.epsilon)
    np.testing.assert_allclose(rh1, rh0, rtol=1e-6)
    # untouched level unchanged in T and q
    np.testing.assert_allclose(out["Q"][0], st["Q"][0], rtol=1e-9)
    # TEOUT = energy of the new state; stratiform state shifted
    te = total_energy(Tn, out["Q"], out["DELP"], st["U"], st["V"], st["PHIS"],
                      out["CLDLIQ"], C)
    np.testing.assert_allclose(out["TEOUT"], te, rtol=1e-12)
    np.testing.assert_allclose(out["TCWAT"] - st["TCWAT"], Tn - T, atol=1e-9)
    np.testing.assert_allclose(out["QCWAT"] - st["QCWAT"], out["Q"] - st["Q"], atol=1e-12)
    # condensate mass conserved
    np.testing.assert_allclose(out["CLDLIQ"] * out["DELP"], st["CLDLIQ"] * st["DELP"])
    assert rep["vapor_added_kg_m2"] > 0 and rep["dPS_mean_Pa"] > 0
    assert rep["te_mismatch_max"] < 1e-12


def test_plugin_refuses_wrong_constants_and_big_vapor_shock():
    st, _ = state()
    wrong = AtmConstants(cpair=1004.0, rair=287.0, zvir=0.608, gravit=9.8, ptop=3.26)
    p = CamAtmPlugin(wrong)
    with pytest.raises(ValueError, match="TEOUT"):
        p.enforce_constraints(st, p.apply_delta(st, np.ones(NL)))
    p = CamAtmPlugin(C)
    with pytest.raises(ValueError, match="vapor column"):
        p.enforce_constraints(st, p.apply_delta(st, np.full(NL, 20.0)))
    ok = CamAtmPlugin(C)
    ok.enforce_constraints(st, ok.apply_delta(st, np.full(NL, 1.0)))


def test_jump_state_rejects_bad_input():
    st, _ = state()
    with pytest.raises(ValueError):
        jump_state(st, np.ones(NL + 1), C)
    with pytest.raises(ValueError):
        jump_state(st, np.full(NL, 30.0), C)
    bad = dict(st, Q=st["Q"].copy())
    bad["Q"][0, 0, 0] = np.nan
    with pytest.raises(ValueError):
        jump_state(bad, np.ones(NL), C)


def test_measured_profile_auto_top_and_taper():
    p = np.geomspace(10, 4e5, 10)
    T_old = np.full(10, 300.0)
    dT = np.array([-3, -2, -1, 0.5, 2, 2.5, 2, 1.5, 1.2, 1.0])
    nocap = dict(ceiling_Pa=1.0, taper_bottom_Pa=2.0)      # ceiling out of the way
    prof = measured_profile(T_old, T_old + 2 * dT, 2.0, p, 10.0,
                            ProfileConfig(smooth=0, **nocap))
    np.testing.assert_allclose(prof.raw_gain, dT)
    assert prof.top_index == 3
    # taper: the top level's gain ramped to zero over two levels
    np.testing.assert_allclose(prof.gain[:3], [0.0, 0.5 / 3, 0.5 * 2 / 3])
    np.testing.assert_allclose(prof.gain[3:], dT[3:])
    fixed = measured_profile(T_old, T_old + 2 * dT, 2.0, p, 10.0,
                             ProfileConfig(smooth=0, top=1e4, taper_levels=0, **nocap))
    assert fixed.gain[p < 1e4].max() == 0 and fixed.gain[-1] == pytest.approx(1.0)
    with pytest.raises(ValueError, match="noise"):
        measured_profile(T_old, T_old + 0.1, 0.1, p, 10.0)


def test_profile_troposphere_ceiling():
    # warming measured all the way up (no auto top): the 100 hPa ceiling cuts it
    p = np.geomspace(100, 4e5, 30)
    T_old = np.full(30, 300.0)
    prof = measured_profile(T_old, T_old + 2.0, 1.0, p, 10.0, ProfileConfig(smooth=0))
    assert np.all(prof.gain[p <= 1e4] == 0.0)                 # nothing at or above 100 hPa
    np.testing.assert_allclose(prof.gain[p >= 2e4], 2.0)      # full strength below 200 hPa
    mid = (p > 1e4) & (p < 2e4)
    assert mid.any() and np.all((prof.gain[mid] > 0) & (prof.gain[mid] < 2.0))
    assert np.all(np.diff(prof.gain[mid]) > 0)                # ramps up with pressure
    assert prof.p_mid[prof.top_index] >= 2e4
    assert prof.summary()["gain_above_ceiling_max"] == 0.0
    with pytest.raises(ValueError):
        ProfileConfig(ceiling_Pa=2e4, taper_bottom_Pa=1e4)
