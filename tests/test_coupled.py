"""Coupled ocean + atmosphere jump: atm-profile -> jump -> check/rollback."""

import json

import numpy as np
import pytest
from synth import gregory_columns, truncate

netCDF4 = pytest.importorskip("netCDF4")

from exocam_accelerate import restart, runstate  # noqa: E402
from exocam_accelerate.atmos import STATE_FIELDS, temperature  # noqa: E402
from exocam_accelerate.cli import main  # noqa: E402
from exocam_accelerate.ocean_advise import advise_ocean  # noqa: E402
from test_atmos import C, NI, NJ, NL, state  # noqa: E402

CASE = "hot"
D0, D_OLD = "0031-01-01-00000", "0021-01-01-00000"
LOG = (f" CPDAIR:     {C.cpair}\n RAIR:       {C.rair}\n ZVIR:      {C.zvir}\n"
       f" SURFACE GRAVITY (m/s2):    {C.gravit}\n Using hyai & hybi from IC:KS= 32  "
       f"PTOP=   {C.ptop}\n")


def write_cam_r(path, st):
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("lev", NL)
        ds.createDimension("pbuf_00008", NL)
        ds.createDimension("lat", NJ)
        ds.createDimension("lon", NI)
        for k, v in st.items():
            dims = (("lat", "lon") if v.ndim == 2 else
                    (("pbuf_00008" if k in ("TCWAT", "QCWAT", "LCWAT", "T_TTEND") else "lev"),
                     "lat", "lon"))
            x = ds.createVariable(k, "f8", dims)
            x[:] = v


def write_trend_series(tdir, T_mean, ps_mean, years=30, dTS_per_yr=0.1):
    """exocam-trend output through model year ``years``: TS in _cam.txt and the
    per-level T / PMID series, warming dTS_per_yr at the surface with a gain of
    0.2 (top) .. 1.0 (bottom) per K, plus a seasonal cycle and weather."""
    tdir.mkdir()
    span = f"{CASE}_0001-01-{years:04d}-12"
    m = np.arange(1, 12 * years + 1)
    t = (m - 0.5) / 12.0
    rng = np.random.default_rng(3)
    season = np.sin(2 * np.pi * t)
    ts = 340.0 + dTS_per_yr * t + season + 0.05 * rng.standard_normal(m.size)
    with open(tdir / f"{span}_cam.txt", "w") as f:
        print("month  TS_native  TS_int1  TS_int2", file=f)
        for i, v in zip(m, ts):
            print(f"{i}  {v:.6f}  {v:.6f}  {v:.6f}", file=f)
    gain = np.linspace(0.2, 1.0, NL)
    T = (T_mean[None, :] + gain[None, :] * dTS_per_yr * (t[:, None] - years)
         + season[:, None] + 0.05 * rng.standard_normal((m.size, NL)))
    P = np.broadcast_to(np.linspace(0.02, 1, NL) * ps_mean, (m.size, NL))
    for name, arr in (("T", T), ("PMID", P)):
        with open(tdir / f"{span}_camlev_{name}.txt", "w") as f:
            print("month  " + "  ".join(f"L{k + 1:02d}" for k in range(NL)), file=f)
            for i, row in zip(m, arr):
                print(f"{i}  " + "  ".join(f"{v:.7g}" for v in row), file=f)


def write_docn_r(path, val):
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("gsize", NJ * NI)
        v = ds.createVariable("somtp", "f8", ("gsize",))
        v[:] = np.full(NJ * NI, val)


@pytest.fixture
def coupled(tmp_path):
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
    write_trend_series(tmp_path / "trend", T.mean((1, 2)), float(st["PS"].mean()))
    return run, arch, st, T


def test_coupled_cli_cycle(coupled, tmp_path, capsys):
    run, arch, st, T = coupled
    adv = advise_ocean(truncate(gregory_columns(), 30), CASE).to_dict()
    af, cf = tmp_path / "a.json", tmp_path / "c.json"
    af.write_text(json.dumps(adv))
    assert main(["atm-profile", "--advice", str(af), "--trend-dir", str(tmp_path / "trend"),
                 "--json", str(cf)]) == 0
    c = json.loads(cf.read_text())
    assert c["schema_version"] == "ocean-atm-1" and c["config"]["heat_ratio"] == 1.0
    assert c["atmosphere"]["method"] == "trend" and c["atmosphere"]["years"] == [21, 30]
    assert len(c["atmosphere"]["sources"]) == 2
    dTS = c["somtp_dT"]
    assert dTS == pytest.approx(adv["dTS_target"])
    assert main(["jump", "--rundir", str(run), "--advice", str(cf), "--archive", str(arch),
                 "--skip-slurm-check", "--yes"]) == 0
    out = capsys.readouterr().out
    assert "cam.r written and verified" in out, out
    cam = run / f"{CASE}.cam.r.{D0}.nc"
    a = restart.read_cam_state(cam)
    Tn = temperature(a["PT"], a["Q"], a["DELP"], C)
    # measured gain 0.2..1.0 per K (profile built from -2 K at the surface)
    np.testing.assert_allclose((Tn - T)[-1], dTS * 1.0, rtol=0.05)
    assert np.allclose(restart.read_somtp(run / f"{CASE}.docn.r.{D0}.nc"), 345.0 + dTS)
    grp = runstate.jump_log_group(run, None)
    assert [json.loads(p.read_text())["plugin"] for p in grp] == ["som_ocean", "cam_atm"]
    assert main(["rollback", "--rundir", str(run), "--archive", str(arch),
                 "--skip-slurm-check", "--yes"]) == 0
    b = restart.read_cam_state(cam)
    assert all(np.array_equal(b[k], st[k]) for k in ("PT", "Q", "DELP", "TEOUT"))
    assert np.allclose(restart.read_somtp(run / f"{CASE}.docn.r.{D0}.nc"), 345.0)
    assert not restart.find_jump_logs(run)


def test_coupled_refuses_mismatched_profile_date(coupled, tmp_path):
    run, arch, _, _ = coupled
    adv = advise_ocean(truncate(gregory_columns(), 30), CASE).to_dict()
    af, cf = tmp_path / "a.json", tmp_path / "c.json"
    af.write_text(json.dumps(adv))
    main(["atm-profile", "--advice", str(af), "--trend-dir", str(tmp_path / "trend"),
          "--json", str(cf)])
    c = json.loads(cf.read_text())
    c["atmosphere"]["restart_date"] = D_OLD
    cf.write_text(json.dumps(c))
    rc = main(["jump", "--rundir", str(run), "--advice", str(cf), "--archive", str(arch),
               "--skip-slurm-check", "--yes"])
    assert rc == 3
    assert not restart.is_jumped(run / f"{CASE}.cam.r.{D0}.nc")


def test_constants_from_log(coupled):
    run, *_ = coupled
    c, f = restart.atm_constants_from_log(run)
    assert c.cpair == pytest.approx(C.cpair) and c.ptop == pytest.approx(C.ptop)
