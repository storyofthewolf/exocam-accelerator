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


def write_cam_i(path, T, ps):
    with netCDF4.Dataset(path, "w") as ds:
        for d, n in (("time", 1), ("lev", NL), ("lat", NJ), ("lon", NI)):
            ds.createDimension(d, n)
        v = ds.createVariable("T", "f8", ("time", "lev", "lat", "lon"))
        v[:] = T[None]
        v = ds.createVariable("PS", "f8", ("time", "lat", "lon"))
        v[:] = ps[None]
        v = ds.createVariable("gw", "f8", ("lat",))
        v[:] = np.ones(NJ)
        for name, arr in (("hyam", np.zeros(NL)), ("hybm", np.linspace(0.02, 1, NL))):
            v = ds.createVariable(name, "f8", ("lev",))
            v[:] = arr
        v = ds.createVariable("P0", "f8", ())
        v[...] = 4e5


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
    for date, dT, somtp in ((D0, 0.0, 345.0), (D_OLD, -2.0, 343.0)):
        rest = arch / "rest" / date
        rest.mkdir(parents=True)
        write_cam_i(rest / f"{CASE}.cam.i.{date}.nc",
                    T + dT * np.linspace(0.2, 1.0, NL)[:, None, None], st["PS"])
        write_docn_r(rest / f"{CASE}.docn.r.{date}.nc", somtp)
    write_cam_r(arch / "rest" / D0 / cam, st)
    return run, arch, st, T


def test_coupled_cli_cycle(coupled, tmp_path, capsys):
    run, arch, st, T = coupled
    adv = advise_ocean(truncate(gregory_columns(), 30), CASE).to_dict()
    af, cf = tmp_path / "a.json", tmp_path / "c.json"
    af.write_text(json.dumps(adv))
    assert main(["atm-profile", "--advice", str(af), "--archive", str(arch),
                 "--json", str(cf)]) == 0
    c = json.loads(cf.read_text())
    assert c["schema_version"] == "ocean-atm-1" and c["config"]["heat_ratio"] == 1.0
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
    main(["atm-profile", "--advice", str(af), "--archive", str(arch), "--json", str(cf)])
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
