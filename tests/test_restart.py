import json

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")

from exocam_accelerate import restart  # noqa: E402
from exocam_accelerate.advise import ADVICE_SCHEMA_VERSION  # noqa: E402
from exocam_accelerate.cli import main  # noqa: E402

CASE = "case"
DATE = "0101-01-01-00000"
NAME = f"{CASE}.cice.r.{DATE}.nc"


def make_cice_r(path, seed=0):
    rng = np.random.default_rng(seed)
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("ni", 6)
        ds.createDimension("nj", 4)
        ds.createDimension("ncat", 5)
        ds.createDimension("ntilyr", 20)
        ds.createDimension("ntslyr", 5)
        for name, dim, lo, hi in [("aicen", "ncat", 0, 0.2), ("vicen", "ncat", 0, 30),
                                  ("vsnon", "ncat", 0, 0.2), ("eicen", "ntilyr", -1e9, -1e8),
                                  ("esnon", "ntslyr", -1e7, -1e6)]:
            v = ds.createVariable(name, "f8", (dim, "nj", "ni"))
            v[:] = rng.uniform(lo, hi, (len(ds.dimensions[dim]), 4, 6))
        ds.istep1 = 87600
    return path


def read(path, name):
    with netCDF4.Dataset(path) as ds:
        return np.array(ds.variables[name][:])


@pytest.fixture
def rundir(tmp_path):
    run = tmp_path / "run"
    run.mkdir()
    make_cice_r(run / NAME)
    (run / "rpointer.ice").write_text(f"./{NAME}\n")
    return run


def test_locate_via_rpointer(rundir):
    assert restart.locate_cice_restart(rundir) == (rundir / NAME).resolve()


def test_bookkeeping_lives_outside_st_archive_globs(rundir):
    # st_archive sweeps ${CASE}.cice.r.[0-9]* in the run dir; ours must not match
    b, log = restart.backup_path(rundir / NAME), restart.log_path(rundir / NAME)
    assert b.parent.name == restart.STATE_DIR
    assert log.parent.name == restart.STATE_DIR
    assert b.name == f"{CASE}.cice.r.{DATE}.pre-accel.nc"


def test_dates_and_years():
    assert restart.restart_date(NAME) == DATE
    assert restart.case_of(NAME) == CASE
    assert restart.first_model_year("0101-01-01-00000") == 101
    assert restart.first_model_year("0101-07-01-00000") == 102


def test_jump_scales_ice_and_keeps_area(rundir):
    path = rundir / NAME
    v0, e0, a0, s0 = (read(path, n) for n in ("vicen", "eicen", "aicen", "vsnon"))
    rec = restart.apply_ice_jump(path, 1.5, provenance={"case": CASE})
    assert rec.written
    np.testing.assert_allclose(read(path, "vicen"), 1.5 * v0)
    np.testing.assert_allclose(read(path, "eicen"), 1.5 * e0)
    np.testing.assert_array_equal(read(path, "aicen"), a0)
    np.testing.assert_array_equal(read(path, "vsnon"), s0)
    np.testing.assert_array_equal(read(rec.backup, "vicen"), v0)
    with netCDF4.Dataset(path) as ds:
        meta = json.loads(ds.getncattr(restart.ATTR))
        assert ds.istep1 == 87600
    assert meta["ice_factor"] == 1.5
    assert meta["jump_model_year"] == 101
    log = json.loads(rec.log.read_text())
    assert log["advice"] == {"case": CASE}
    assert log["restart_date"] == DATE
    assert restart.find_jump_logs(rundir) == [rec.log]


def test_rejump_does_not_compound(rundir):
    path = rundir / NAME
    v0 = read(path, "vicen")
    restart.apply_ice_jump(path, 1.5)
    restart.apply_ice_jump(path, 1.2)
    np.testing.assert_allclose(read(path, "vicen"), 1.2 * v0)


def test_dry_run_writes_nothing(rundir):
    path = rundir / NAME
    v0 = read(path, "vicen")
    rec = restart.apply_ice_jump(path, 1.5, dry_run=True)
    assert not rec.written
    assert not rec.backup.exists()
    np.testing.assert_array_equal(read(path, "vicen"), v0)


def test_refuses_jumped_file_without_backup(rundir):
    path = rundir / NAME
    rec = restart.apply_ice_jump(path, 1.5)
    rec.backup.unlink()
    with pytest.raises(RuntimeError, match="scale twice"):
        restart.apply_ice_jump(path, 1.5)


def test_restore_retires_log(rundir):
    path = rundir / NAME
    v0 = read(path, "vicen")
    rec = restart.apply_ice_jump(path, 1.5)
    restart.restore(path)
    np.testing.assert_array_equal(read(path, "vicen"), v0)
    assert not rec.log.exists()
    assert rec.log.with_name(rec.log.name[:-5] + ".restored.json").exists()
    assert restart.find_jump_logs(rundir) == []


def test_snow_factor(rundir):
    path = rundir / NAME
    s0 = read(path, "vsnon")
    restart.apply_ice_jump(path, 1.5, snow_factor=1.1)
    np.testing.assert_allclose(read(path, "vsnon"), 1.1 * s0)


def test_cli_jump_from_advice(rundir, tmp_path):
    path = rundir / NAME
    v0 = read(path, "vicen")
    advice = tmp_path / "advice.json"
    advice.write_text(json.dumps({"case": CASE, "ice_factor": 1.3, "model_year": 100,
                                  "schema_version": ADVICE_SCHEMA_VERSION}))
    assert main(["jump", "--rundir", str(rundir), "--advice", str(advice),
                 "--allow-no-archive-rollback", "--skip-slurm-check", "--yes"]) == 0
    np.testing.assert_allclose(read(path, "vicen"), 1.3 * v0)
    assert main(["restore", "--rundir", str(rundir)]) == 0
    np.testing.assert_array_equal(read(path, "vicen"), v0)


def test_cli_jump_blocked_without_squeue(rundir, monkeypatch):
    import exocam_accelerate.runstate as rs
    monkeypatch.setattr(rs, "active_jobs", lambda: None)
    v0 = read(rundir / NAME, "vicen")
    assert main(["jump", "--rundir", str(rundir), "--ice-factor", "1.2", "--yes"]) == 3
    np.testing.assert_array_equal(read(rundir / NAME, "vicen"), v0)


def test_cli_refuses_no_jump_advice(rundir, tmp_path):
    advice = tmp_path / "advice.json"
    advice.write_text(json.dumps({"case": CASE, "ice_factor": None,
                                  "reasons": ["ice edge still moving"]}))
    assert main(["jump", "--rundir", str(rundir), "--advice", str(advice),
                 "--skip-slurm-check", "--yes"]) == 1


def test_cli_hard_bound(rundir):
    assert main(["jump", "--rundir", str(rundir), "--ice-factor", "3",
                 "--allow-no-archive-rollback", "--skip-slurm-check", "--yes"]) == 2


# ---- tapered jump ----

OLD_DATE = "0091-01-01-00000"


def _set(path, **fields):
    with netCDF4.Dataset(path, "r+") as ds:
        for name, fn in fields.items():
            ds.variables[name][:] = fn(np.array(ds.variables[name][:]))


@pytest.fixture
def archive(rundir, tmp_path):
    """Pristine archived restarts at DATE and 10 years earlier, plus a grid file."""
    import shutil
    _set(rundir / NAME, aicen=lambda a: np.full_like(a, 0.2))
    arch = tmp_path / "archive"
    now = arch / "rest" / DATE
    old = arch / "rest" / OLD_DATE
    now.mkdir(parents=True)
    old.mkdir(parents=True)
    shutil.copy2(rundir / NAME, now / NAME)
    old_name = f"{CASE}.cice.r.{OLD_DATE}.nc"
    shutil.copy2(rundir / NAME, old / old_name)
    # ice grew by 25 % everywhere except column 0, which sat still
    # (near local equilibrium)
    _set(old / old_name, vicen=lambda v: np.concatenate(
        [v[:, :, :1], 0.8 * v[:, :, 1:]], axis=2))
    hist = arch / "ice" / "hist"
    hist.mkdir(parents=True)
    with netCDF4.Dataset(hist / f"{CASE}.cice.h.0100-12.nc", "w") as ds:
        ds.createDimension("nj", 4)
        ds.createDimension("ni", 6)
        ds.createVariable("tarea", "f4", ("nj", "ni"))[:] = np.ones((4, 6))
        ds.createVariable("tmask", "f4", ("nj", "ni"))[:] = np.ones((4, 6))
    return arch


def _uniform_advice(tmp_path):
    from synth import stefan_columns, truncate
    from exocam_accelerate.advise import AdvisorConfig, advise
    adv = advise(truncate(stefan_columns(), 60), CASE,
                 AdvisorConfig(max_ice_factor=1.5)).to_dict()
    adv["model_year"] = 100
    p = tmp_path / "advice.json"
    p.write_text(json.dumps(adv))
    return p


def test_build_taper_mask_from_archive(archive):
    mask, sources = restart.build_taper_mask(archive, CASE, DATE, 10)
    assert mask.weight.shape == (4, 6)
    assert np.all(mask.weight[:, 0] == 0.0)          # the static column
    assert mask.weight[:, 1:].max() == 1.0
    assert len(sources) == 3


def test_taper_mask_refuses_jumped_baseline(archive):
    restart.apply_ice_jump(archive / "rest" / DATE / NAME, 1.2)
    with pytest.raises(RuntimeError, match="pristine"):
        restart.build_taper_mask(archive, CASE, DATE, 10)


def test_cli_taper_then_jump(rundir, archive, tmp_path):
    from exocam_accelerate.advise import TAPERED_ADVICE_SCHEMA_VERSION
    from exocam_accelerate.taper import cell_factors
    path = rundir / NAME
    v0, e0 = read(path, "vicen"), read(path, "eicen")
    out = tmp_path / "tapered.json"
    assert main(["taper", "--advice", str(_uniform_advice(tmp_path)),
                 "--archive", str(archive), "--json", str(out)]) == 0
    adv = json.loads(out.read_text())
    assert adv["schema_version"] == TAPERED_ADVICE_SCHEMA_VERSION
    weights = tmp_path / "tapered.taper.nc"
    assert adv["taper"]["weights_file"] == str(weights.resolve())
    w = restart.read_taper_file(weights)["weight"]

    assert main(["jump", "--rundir", str(rundir), "--advice", str(out),
                 "--allow-no-archive-rollback", "--skip-slurm-check", "--yes"]) == 0
    f = cell_factors(adv["ice_factor"], w)
    np.testing.assert_allclose(read(path, "vicen"), v0 * f)
    np.testing.assert_allclose(read(path, "eicen"), e0 * f)
    np.testing.assert_allclose(read(path, "vicen")[:, :, 0], v0[:, :, 0])
    log = json.loads(restart.log_path(path).read_text())
    assert log["taper"]["sha256"] == adv["taper"]["weights_sha256"]
    assert (restart.state_dir(path) / (NAME + restart.TAPER_SUFFIX)).is_file()


def test_cli_jump_refuses_altered_weight_map(rundir, archive, tmp_path):
    out = tmp_path / "tapered.json"
    assert main(["taper", "--advice", str(_uniform_advice(tmp_path)),
                 "--archive", str(archive), "--json", str(out)]) == 0
    _set(tmp_path / "tapered.taper.nc", weight=lambda w: np.ones_like(w))
    v0 = read(rundir / NAME, "vicen")
    assert main(["jump", "--rundir", str(rundir), "--advice", str(out),
                 "--allow-no-archive-rollback", "--skip-slurm-check", "--yes"]) == 2
    np.testing.assert_array_equal(read(rundir / NAME, "vicen"), v0)
