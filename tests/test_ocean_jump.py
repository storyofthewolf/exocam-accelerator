"""docn.r somtp jump: restart layer, post-jump check, CLI, pre-flight, rollback."""

import json

import numpy as np
import pytest
from synth import gregory_columns, truncate

netCDF4 = pytest.importorskip("netCDF4")

from exocam_accelerate import restart, runstate  # noqa: E402
from exocam_accelerate.check import Verdict, check_any, check_ocean_jump  # noqa: E402
from exocam_accelerate.cli import main  # noqa: E402
from exocam_accelerate.ocean_advise import advise_ocean  # noqa: E402
from exocam_accelerate.som_ocean import TK_FRZ_SW  # noqa: E402

CASE = "hot"
D0 = "0031-01-01-00000"
D_OLD = "0021-01-01-00000"
NJ, NI = 4, 6
JUMP_YEAR = 31


def make_docn_r(path, somtp):
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("gsize", somtp.size)
        v = ds.createVariable("somtp", "f8", ("gsize",))
        v[:] = somtp
        ds.file_version = "shr_pcdf_v0_0_01"
    return path


def make_domain(path, land=None):
    lat = np.linspace(-60, 60, NJ)
    lon = np.arange(NI) * 60.0
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("nj", NJ)
        ds.createDimension("ni", NI)
        mask = np.ones((NJ, NI))
        if land is not None:
            mask[land] = 0
        for name, arr in (("mask", mask),
                          ("area", np.cos(np.radians(lat))[:, None] * np.ones((1, NI))),
                          ("xc", lon[None, :] * np.ones((NJ, 1))),
                          ("yc", lat[:, None] * np.ones((1, NI)))):
            v = ds.createVariable(name, "f8", ("nj", "ni"))
            v[:] = arr
    return path


def somtp_field(base=340.0, grad=10.0):
    lat = np.linspace(-60, 60, NJ)
    return (base + grad * np.cos(np.radians(lat))[:, None] * np.ones((1, NI))).ravel()


@pytest.fixture
def case(tmp_path):
    """A run directory at D0 with an archived set at D0 and D_OLD."""
    run, arch = tmp_path / "run", tmp_path / "archive"
    run.mkdir()
    name = f"{CASE}.docn.r.{D0}.nc"
    make_docn_r(run / name, somtp_field())
    (run / "rpointer.ocn").write_text(f"{name}\n{CASE}.docn.rs1.{D0}.bin\n")
    (run / "rpointer.atm").write_text(f"{CASE}.cam.r.{D0}.nc\n")
    for date, field in ((D0, somtp_field()), (D_OLD, somtp_field(335.0, 5.0))):
        rest = arch / "rest" / date
        rest.mkdir(parents=True)
        make_docn_r(rest / f"{CASE}.docn.r.{date}.nc", field)
    dom = make_domain(tmp_path / "domain.nc")
    return run, arch, run / name, dom


def advice_dict():
    a = advise_ocean(truncate(gregory_columns(), JUMP_YEAR - 1), CASE)
    assert a.jump, a.reasons
    return a.to_dict()


# ---------------------------------------------------------------- restart
def test_locate_and_case(case):
    run, _, path, _ = case
    assert restart.locate_docn_restart(run) == path.resolve()
    assert restart.case_of(path) == CASE


def test_jump_backup_log_and_reapply(case):
    run, _, path, dom = case
    before = restart.read_somtp(path)
    rec = restart.apply_ocean_jump(path, 2.5, {"note": 1}, domain_file=dom)
    np.testing.assert_allclose(restart.read_somtp(path), before + 2.5)
    assert rec.backup.parent.name == restart.STATE_DIR
    assert restart.is_jumped(path) and not restart.is_jumped(rec.backup)
    log = json.loads(rec.log.read_text())
    assert log["plugin"] == "som_ocean" and log["somtp_dT_applied_mean"] == pytest.approx(2.5)
    assert log["jump_model_year"] == JUMP_YEAR
    # re-applying reads the pristine copy: no compounding
    restart.apply_ocean_jump(path, 1.0, domain_file=dom)
    np.testing.assert_allclose(restart.read_somtp(path), before + 1.0)
    restart.restore(path)
    np.testing.assert_allclose(restart.read_somtp(path), before)


def test_masked_and_frozen_cells_untouched(case, tmp_path):
    run, _, path, _ = case
    t = somtp_field()
    t[0] = TK_FRZ_SW
    make_docn_r(path, t)
    dom = make_domain(tmp_path / "d2.nc", land=(1, 2))
    rec = restart.apply_ocean_jump(path, 3.0, domain_file=dom)
    after = restart.read_somtp(path)
    assert after[0] == TK_FRZ_SW
    assert after[1 * NI + 2] == t[1 * NI + 2]
    assert rec.metadata["somtp_dT_applied_mean"] == pytest.approx(3.0)


def test_refuses_double_jump_without_backup(case):
    _, _, path, _ = case
    restart.apply_ocean_jump(path, 1.0)
    restart.backup_path(path).unlink()
    with pytest.raises(RuntimeError, match="scale|shift"):
        restart.apply_ocean_jump(path, 1.0)


def test_dry_run_writes_nothing(case):
    _, _, path, _ = case
    before = restart.read_somtp(path)
    rec = restart.apply_ocean_jump(path, 2.0, dry_run=True)
    assert not rec.written and not rec.backup.exists()
    np.testing.assert_array_equal(restart.read_somtp(path), before)


def test_pattern_and_map(case, tmp_path):
    run, arch, path, dom = case
    pat, grid, sources = restart.build_ocean_pattern(arch, CASE, D0, dom, 10)
    f = tmp_path / "p.pattern.nc"
    restart.write_pattern_file(f, pat, grid, 2.0)
    w = restart.read_pattern_file(f)["weight"]
    assert w.shape == (NJ, NI)
    before = restart.read_somtp(path)
    rec = restart.apply_ocean_jump(path, 2.0, pattern_file=f, domain_file=dom)
    d = restart.read_somtp(path) - before
    np.testing.assert_allclose(d, 2.0 * w.ravel())
    assert rec.metadata["somtp_dT_applied_mean"] == pytest.approx(2.0)
    out = restart.write_somtp_map(path, dom, tmp_path / "m.nc")
    with netCDF4.Dataset(out) as ds:
        assert ds.variables["somtp"].dimensions == ("lat", "lon")
        np.testing.assert_allclose(ds.variables["somtp"][:], (before + d).reshape(NJ, NI))


def test_domain_size_mismatch(case, tmp_path):
    _, _, path, _ = case
    make_docn_r(path, np.full(7, 300.0))
    with pytest.raises(ValueError, match="cells"):
        restart.apply_ocean_jump(path, 1.0, domain_file=make_domain(tmp_path / "x.nc"))


# ---------------------------------------------------------------- check
def jump_log(dT=None):
    adv = advice_dict()
    dT = adv["somtp_dT"] if dT is None else dT
    return {"plugin": "som_ocean", "jump_model_year": JUMP_YEAR, "somtp_dT": dT,
            "somtp_dT_applied_mean": dT, "advice": adv}


def run_after(years_after, dT=None, **kw):
    log = jump_log(dT)
    cols = gregory_columns(years=JUMP_YEAR - 1 + years_after, jump_year=JUMP_YEAR,
                           dT=log["somtp_dT"], **kw)
    return check_ocean_jump(cols, log), log


def test_pass_when_run_follows_the_line():
    res, _ = run_after(6)
    assert res.verdict is Verdict.PASS, res.reasons
    assert abs(res.metrics["dN"]) < 0.05
    assert res.metrics["land_fraction"] == pytest.approx(1.0, abs=0.35)
    assert check_any(truncate(gregory_columns(years=36, jump_year=JUMP_YEAR,
                                              dT=res.metrics["somtp_dT_applied"]), 36),
                     jump_log()).verdict is Verdict.PASS


def test_wait_then_settling():
    log = jump_log()
    assert check_ocean_jump(truncate(gregory_columns(), JUMP_YEAR - 1),
                            log).verdict is Verdict.WAIT
    res, _ = run_after(1)
    assert res.verdict is Verdict.WAIT


def test_fail_when_jump_never_landed():
    log = jump_log()
    cols = gregory_columns(years=JUMP_YEAR + 4)          # no jump in the data
    res = check_ocean_jump(cols, log)
    assert res.verdict is Verdict.FAIL
    assert any("did not land" in r for r in res.reasons)


def test_fail_when_off_the_line():
    res, _ = run_after(6, N_offset_after=3.0)
    assert res.verdict is Verdict.FAIL
    assert any("off the Gregory line" in r for r in res.reasons)


def test_fail_on_nan_after_jump():
    log = jump_log()
    cols = gregory_columns(years=JUMP_YEAR + 5, jump_year=JUMP_YEAR, dT=log["somtp_dT"])
    cols["TS_native"] = cols["TS_native"].copy()
    cols["TS_native"][-12:] = np.nan
    assert check_ocean_jump(cols, log).verdict is Verdict.FAIL


def test_heat_ratio_implied_measures_the_shortfall():
    # the model kept only 2/3 of the increment (the atmosphere took the rest)
    log = jump_log()
    cols = gregory_columns(years=JUMP_YEAR + 5, jump_year=JUMP_YEAR,
                           dT=log["somtp_dT"] * 2 / 3)
    res = check_ocean_jump(cols, log)
    assert res.metrics["heat_ratio_implied"] == pytest.approx(1.5, rel=0.15)
    assert res.verdict is Verdict.PASS          # still on the line: not a failure


# ---------------------------------------------------------------- CLI + pre-flight
def test_cli_jump_check_rollback(case, tmp_path, monkeypatch, capsys):
    run, arch, path, dom = case
    adv = advice_dict()
    af = tmp_path / "adv.json"
    af.write_text(json.dumps(adv))
    before = restart.read_somtp(path)
    rc = main(["jump", "--rundir", str(run), "--advice", str(af), "--archive", str(arch),
               "--domain-file", str(dom), "--skip-slurm-check", "--yes"])
    out = capsys.readouterr().out
    assert rc == 0, out
    np.testing.assert_allclose(restart.read_somtp(path), before + adv["somtp_dT"])
    assert "rpointer" in out

    rc = main(["restore", "--rundir", str(run), "--ocean"])
    assert rc == 0
    np.testing.assert_allclose(restart.read_somtp(path), before)

    # rollback from the archive restores the pristine docn.r
    main(["jump", "--docn-r", str(path), "--delta-t", "1.5", "--archive", str(arch),
          "--skip-slurm-check", "--yes"])
    rc = main(["rollback", "--rundir", str(run), "--archive", str(arch),
               "--skip-slurm-check", "--yes"])
    assert rc == 0, capsys.readouterr().out
    np.testing.assert_allclose(restart.read_somtp(path), before)


def test_preflight_blocks_mismatched_rpointer(case):
    run, arch, path, _ = case
    (run / "rpointer.atm").write_text(f"{CASE}.cam.r.0021-01-01-00000.nc\n")
    f = runstate.preflight(path, arch, None, advice_dict())
    assert runstate.blocked(f)


def test_preflight_accepts_ocean_schema(case):
    run, arch, path, _ = case
    f = runstate.preflight(path, arch, None, advice_dict())
    assert not runstate.blocked(f), [x.message for x in f if x.level == "block"]


def test_cli_advise_ocean_and_map(tmp_path, capsys):
    cols = gregory_columns(years=30)
    hdr = ["month"] + [k for k in cols if k != "month"]
    rows = np.column_stack([cols[k] for k in hdr])
    p = tmp_path / f"{CASE}_0001-01-0030-12_cam.txt"
    np.savetxt(p, rows, header="  ".join(hdr), comments="", fmt="%.6f")
    rc = main(["advise-ocean", str(tmp_path), CASE, "--json", str(tmp_path / "a.json")])
    out = capsys.readouterr().out
    assert rc == 0 and "RECOMMENDED somtp increment" in out, out
    assert json.loads((tmp_path / "a.json").read_text())["plugin"] == "som_ocean"


def test_cli_somtp_map(case, tmp_path, capsys):
    run, _, path, dom = case
    rc = main(["somtp-map", "--rundir", str(run), "--domain-file", str(dom),
               "-o", str(tmp_path / "m.nc")])
    assert rc == 0 and (tmp_path / "m.nc").is_file()


def test_domain_from_docn_ocn_in(case, tmp_path, capsys):
    run, _, path, dom = case
    (run / "docn_ocn_in").write_text(f"&shr_strdata_nml\n  datamode = 'SOM'\n"
                                     f"  domainfile = '{dom}'\n/\n")
    assert restart.find_docn_domain(run) == dom
    rc = main(["somtp-map", "--rundir", str(run), "-o", str(tmp_path / "m2.nc")])
    assert rc == 0


# ---------------------------------------------------------------- probe check
from exocam_accelerate.check import check_ocean_probe  # noqa: E402
from exocam_accelerate.ocean_advise import ProbeConfig, probe_ocean  # noqa: E402


def probe_log(dT, tau=40.0):
    adv = probe_ocean(truncate(gregory_columns(tau=tau, noise=0.02), JUMP_YEAR - 1),
                      CASE, ProbeConfig(probe_dT=dT))
    assert adv["somtp_dT"] is not None, adv["reasons"]
    return {"plugin": "som_ocean", "jump_model_year": JUMP_YEAR, "somtp_dT": dT,
            "somtp_dT_applied_mean": dT, "advice": adv}


def probe_run(dT, years_after, tau=40.0, **kw):
    log = probe_log(dT, tau)
    cols = gregory_columns(years=JUMP_YEAR - 1 + years_after, tau=tau, jump_year=JUMP_YEAR,
                           dT=dT, noise=0.02, **kw)
    return check_any(cols, log), log


def test_probe_reads_lambda_and_side():
    res, _ = probe_run(5.0, 8)
    assert res.verdict is Verdict.PASS, res.reasons
    assert res.metrics["lambda"] == pytest.approx(-1 / -1.2 * 1.0, rel=0.15)   # -1/c1
    assert res.metrics["TS_eq"] == pytest.approx(350.0, abs=1.0)
    assert res.metrics["bracketed"] == 0.0                   # still below equilibrium


def test_probe_overshoot_is_bracketed():
    res, _ = probe_run(25.0, 8, tau=40.0)
    assert res.verdict is Verdict.PASS, res.reasons
    assert res.metrics["bracketed"] == 1.0
    assert res.metrics["TS_eq"] == pytest.approx(350.0, abs=1.5)


def test_probe_runaway_signature_fails():
    # the imbalance rises after the warm probe and TS gives no readout to the
    # contrary (pre-probe window too short for it): no restoring feedback
    log = probe_log(5.0)
    log["advice"]["since_year"] = JUMP_YEAR - 4
    cols = gregory_columns(years=JUMP_YEAR - 1 + 8, tau=40.0, jump_year=JUMP_YEAR,
                           dT=5.0, noise=0.02, N_offset_after=12.0)
    res = check_any(cols, log)
    assert res.verdict is Verdict.FAIL
    assert any("runaway" in r for r in res.reasons)


def test_probe_energy_runaway_downgraded_when_ts_restores():
    # same energy signal, but the TS drift clearly shows the probe restoring:
    # an energy false alarm must not trigger a rollback
    res, _ = probe_run(5.0, 8, N_offset_after=12.0)
    assert res.verdict is Verdict.WAIT, res.reasons
    assert any("false alarm" in r for r in res.reasons)


def test_probe_waits_then_not_landed():
    res, _ = probe_run(5.0, 3)
    assert res.verdict is Verdict.WAIT
    log = probe_log(5.0)
    res = check_ocean_probe(gregory_columns(years=JUMP_YEAR + 6, tau=40.0, noise=0.02), log)
    assert res.verdict is Verdict.FAIL and any("did not land" in r for r in res.reasons)


def test_probe_advice_passes_preflight(case):
    run, arch, path, _ = case
    adv = probe_log(3.0)["advice"]
    f = runstate.preflight(path, arch, None, adv)
    assert not runstate.blocked(f), [x.message for x in f if x.level == "block"]
