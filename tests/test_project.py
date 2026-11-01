"""Projection jump: project -> jump --override-limits -> rollback, and the
vapor-scaled atmosphere edit."""

import json

import numpy as np
import pytest

from exocam_accelerate.atmos import (CamAtmPlugin, interfaces, jump_state, p_mid, qsat,
                                     temperature)
from exocam_accelerate.project import (ProjectConfig, linear_projection, project_state,
                                       superadiabatic_levels)
from exocam_accelerate.som_ocean import check_dT
from test_atmos import C, NL, state

KAPPA, EPS = 0.2858, 0.62
NP = 20
P = np.logspace(3, 5.6, NP)                   # 10 .. 3981 hPa, model order (top first)


def column(frac, years=30, dTS_per_yr=0.4, vapor_per_yr=0.01):
    T = 400.0 * (P / P[-1]) ** (frac * KAPPA)
    y = np.arange(1, years + 1)
    TS = T[-1] + dTS_per_yr * (y - years)
    W = 1000.0 * (1.0 + vapor_per_yr * (y - years))
    return (y, TS, y, W, y, np.tile(T, (years, 1)), np.tile(np.full(NP, 0.01), (years, 1)),
            np.tile(P, (years, 1)))


def test_linear_projection():
    y = np.arange(1, 31)
    s, now, tgt = linear_projection(y, 2.0 * y + 1.0, 30, 130, 20)
    assert s == pytest.approx(2.0) and now == 61.0 and tgt == pytest.approx(261.0)
    with pytest.raises(ValueError, match="complete years"):
        linear_projection(y, y, 30, 130, 40)


def test_project_shifts_troposphere_and_scales_vapor():
    r = project_state(*column(0.5), 30, 130, KAPPA, EPS, ProjectConfig())
    assert r["dTS"] == pytest.approx(40.0)
    assert r["vapor_factor"] == pytest.approx(2.0)
    full = P >= 2.0e4
    np.testing.assert_allclose(r["dT_levels"][full], 40.0)
    np.testing.assert_allclose(r["vapor_scale"][full], 2.0)
    assert np.all(r["dT_levels"][P <= r["ceiling_Pa"]] == 0.0)
    assert np.all(r["vapor_scale"][P <= r["ceiling_Pa"]] == 1.0)
    assert not superadiabatic_levels(r["T_after"], P, KAPPA)


def test_project_raises_ceiling_for_lapse_rate():
    r = project_state(*column(0.5), 30, 130, KAPPA, EPS, ProjectConfig())
    assert r["ceiling_Pa"] < 1.0e4
    assert r["adjustments"] and r["adjustments"][0].startswith("lapse rate")
    # an already very stable column needs no adjustment
    s = project_state(*column(0.1), 30, 130, KAPPA, EPS, ProjectConfig())
    assert s["ceiling_Pa"] == 1.0e4 and not s["adjustments"]


def test_project_refuses_k_table_range_and_bad_target():
    with pytest.raises(ValueError, match="k-table"):
        project_state(*column(0.5, dTS_per_yr=1.0), 30, 130, KAPPA, EPS, ProjectConfig())
    with pytest.raises(ValueError, match="not after"):
        project_state(*column(0.5), 30, 30, KAPPA, EPS, ProjectConfig())


def test_jump_state_vapor_scale_keeps_dry_mass_and_caps_at_saturation():
    st, T = state()
    vs = np.full(NL, 1.5)
    out, rep = jump_state(st, np.full(NL, 5.0), C, vapor_scale=vs)
    np.testing.assert_allclose(out["DELP"] * (1 - out["Q"]), st["DELP"] * (1 - st["Q"]),
                               rtol=1e-12)
    r0 = st["Q"] / (1 - st["Q"])
    r1 = out["Q"] / (1 - out["Q"])
    uncapped = r1 < 1.5 * r0 * (1 - 1e-9)
    np.testing.assert_allclose(r1[~uncapped], 1.5 * r0[~uncapped], rtol=1e-9)
    assert rep["humidity"] == "vapor scale" and rep["saturation_capped_cells"] == 0
    # a huge scale is capped at saturation at the new T
    out, rep = jump_state(st, np.full(NL, 5.0), C, vapor_scale=np.full(NL, 50.0))
    Tn = temperature(out["PT"], out["Q"], out["DELP"], C)
    qs = qsat(Tn, p_mid(interfaces(out["DELP"], C.ptop)), C.epsilon)
    assert rep["saturation_capped_cells"] > 0
    assert np.all(out["Q"] <= np.maximum(qs, st["Q"]) * (1 + 1e-6))


def test_override_bounds():
    st, _ = state()
    with pytest.raises(ValueError, match="level increment"):
        jump_state(st, np.full(NL, 40.0), C)
    out, rep = jump_state(st, np.full(NL, 40.0), C, override=True)
    assert rep["override_limits"] and rep["dT_max"] == pytest.approx(40.0)
    with pytest.raises(ValueError, match="k-table"):
        jump_state(st, np.full(NL, 140.0), C, override=True)
    # the vapor-shock refusal is lifted only with override
    p = CamAtmPlugin(C, vapor_scale=np.full(NL, 3.0))
    with pytest.raises(ValueError, match="vapor column"):
        p.enforce_constraints(st, p.apply_delta(st, np.full(NL, 5.0)))
    p = CamAtmPlugin(C, vapor_scale=np.full(NL, 3.0), override=True)
    p.enforce_constraints(st, p.apply_delta(st, np.full(NL, 5.0)))
    with pytest.raises(ValueError, match="hard bound"):
        check_dT(40.0)
    assert check_dT(40.0, override=True) == 40.0


# --- command line: project -> jump -> rollback --------------------------------

netCDF4 = pytest.importorskip("netCDF4")

from exocam_accelerate import restart  # noqa: E402
from exocam_accelerate.cli import main  # noqa: E402
from test_coupled import CASE, D0, LOG, write_cam_r, write_docn_r  # noqa: E402


def write_projection_trends(tdir, T_mean, Q_mean, P_mean, years=30):
    """TS and TMQ in _cam.txt, per-level T / Q / PMID, through model year
    ``years``: the surface warms 0.3 K/yr (a 100-yr projection crosses the
    25 K bound), the vapor column 0.5 %/yr."""
    tdir.mkdir()
    span = f"{CASE}_0001-01-{years:04d}-12"
    m = np.arange(1, 12 * years + 1)
    t = (m - 0.5) / 12.0
    ts = T_mean[-1] + 0.3 * (t - years)
    tmq = 500.0 * (1 + 0.005 * (t - years))
    with open(tdir / f"{span}_cam.txt", "w") as f:
        print("month  TS_native  TS_int1  TS_int2  TMQ_native  TMQ_int1  TMQ_int2", file=f)
        for i, a, b in zip(m, ts, tmq):
            print(f"{i}  {a:.6f}  {a:.6f}  {a:.6f}  {b:.6f}  {b:.6f}  {b:.6f}", file=f)
    rows = {"T": T_mean[None, :] + 0.3 * (t[:, None] - years),
            "Q": np.broadcast_to(Q_mean, (m.size, NL)),
            "PMID": np.broadcast_to(P_mean, (m.size, NL))}
    for name, arr in rows.items():
        with open(tdir / f"{span}_camlev_{name}.txt", "w") as f:
            print("month  " + "  ".join(f"L{k + 1:02d}" for k in range(NL)), file=f)
            for i, row in zip(m, arr):
                print(f"{i}  " + "  ".join(f"{v:.7g}" for v in row), file=f)


@pytest.fixture
def projected(tmp_path):
    run, arch = tmp_path / "run", tmp_path / "archive"
    run.mkdir()
    st, T = state()
    cam, docn = f"{CASE}.cam.r.{D0}.nc", f"{CASE}.docn.r.{D0}.nc"
    write_cam_r(run / cam, st)
    write_docn_r(run / docn, 345.0)
    (run / "rpointer.atm").write_text(cam + "\n")
    (run / "rpointer.ocn").write_text(f"{docn}\n{CASE}.docn.rs1.{D0}.bin\n")
    (run / "atm.log.260930-000000").write_text(LOG)
    rest = arch / "rest" / D0
    rest.mkdir(parents=True)
    write_docn_r(rest / docn, 345.0)
    write_cam_r(rest / cam, st)
    pm = p_mid(interfaces(st["DELP"], C.ptop)).mean((1, 2))
    write_projection_trends(tmp_path / "trend", T.mean((1, 2)), st["Q"].mean((1, 2)), pm)
    return run, arch, st, T


def test_project_cli_cycle(projected, tmp_path, capsys):
    run, arch, st, T = projected
    cf = tmp_path / "p.json"
    assert main(["project", str(tmp_path / "trend"), CASE, "--target-year", "130",
                 "--atm-log", str(run), "--json", str(cf)]) == 0
    c = json.loads(cf.read_text())
    assert c["mode"] == "project" and c["override_limits"] is True
    assert c["somtp_dT"] == pytest.approx(30.0, abs=0.1)
    assert c["vapor_factor"] == pytest.approx(1.5, rel=1e-3)
    assert c["atmosphere"]["restart_date"] == D0
    jump = ["jump", "--rundir", str(run), "--advice", str(cf), "--archive", str(arch),
            "--skip-slurm-check", "--yes"]
    assert main(jump) == 2
    assert "--override-limits must be given" in capsys.readouterr().err
    assert not restart.is_jumped(run / f"{CASE}.cam.r.{D0}.nc")
    assert main(jump + ["--override-limits"]) == 0
    out = capsys.readouterr().out
    assert "OVERRIDE" in out and "vapor scale" in out, out
    a = restart.read_cam_state(run / f"{CASE}.cam.r.{D0}.nc")
    Tn = temperature(a["PT"], a["Q"], a["DELP"], C)
    np.testing.assert_allclose((Tn - T)[-1], c["somtp_dT"], rtol=1e-6)
    assert np.allclose(restart.read_somtp(run / f"{CASE}.docn.r.{D0}.nc"),
                       345.0 + c["somtp_dT"])
    assert main(["rollback", "--rundir", str(run), "--archive", str(arch),
                 "--skip-slurm-check", "--yes"]) == 0
    b = restart.read_cam_state(run / f"{CASE}.cam.r.{D0}.nc")
    assert all(np.array_equal(b[k], st[k]) for k in ("PT", "Q", "DELP", "TEOUT"))


def test_override_flag_refused_for_ordinary_advice(projected, tmp_path, capsys):
    run, arch, _, _ = projected
    assert main(["jump", "--rundir", str(run), "--delta-t", "2", "--archive", str(arch),
                 "--skip-slurm-check", "--yes", "--override-limits"]) == 2
    assert "only to projection" in capsys.readouterr().err
    assert not restart.is_jumped(run / f"{CASE}.docn.r.{D0}.nc")
