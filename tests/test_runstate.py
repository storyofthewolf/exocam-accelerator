"""Pre-flight and rollback against a fake CESM run directory + archive.

The fake mirrors CESM 1.2.1 st_archive: at the end of each segment the newest
restart set and rpointers go to <archive>/rest/<date>/ and are copied back
into the run directory.
"""

import json
import shutil

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")

from test_restart import make_cice_r, read  # noqa: E402

from exocam_accelerate import restart, runstate  # noqa: E402
from exocam_accelerate.advise import ADVICE_SCHEMA_VERSION  # noqa: E402
from exocam_accelerate.cli import main  # noqa: E402

CASE = "case"
D0 = "0101-01-01-00000"      # jump date
D1 = "0106-01-01-00000"      # end of the post-jump segment


def no_jobs():
    return set()


def write_set(run, date, seed=0):
    """A restart set + rpointers for ``date`` in ``run``."""
    make_cice_r(run / f"{CASE}.cice.r.{date}.nc", seed)
    for comp in ("cam.r", "cam.rs", "cpl.r", "docn.r"):
        (run / f"{CASE}.{comp}.{date}.nc").write_text(comp + date)
    (run / "rpointer.ice").write_text(f"{CASE}.cice.r.{date}.nc\n")
    (run / "rpointer.drv").write_text(f"{CASE}.cpl.r.{date}.nc\n")
    (run / "rpointer.ocn").write_text(f"{CASE}.docn.r.{date}.nc\n")
    (run / "rpointer.atm").write_text(
        f"{CASE}.cam.r.{date}.nc\n\n# comment naming {CASE}.cam.r.{date}.nc\n")


def st_archive(run, archive, date):
    rest = archive / "rest" / date
    rest.mkdir(parents=True)
    for p in run.iterdir():
        if p.is_file() and (p.name.startswith("rpointer.") or date in p.name):
            shutil.copy2(p, rest / p.name)


@pytest.fixture
def case(tmp_path):
    run, archive = tmp_path / "run", tmp_path / "archive"
    run.mkdir()
    write_set(run, D0)
    for m in range(1, 13):
        (run / f"{CASE}.cam.h0.0100-{m:02d}.nc").write_text("h0")
    st_archive(run, archive, D0)
    return run, archive


def levels(findings):
    return {f.level for f in findings}


class TestPreflight:
    def test_clean_case_passes(self, case):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, no_jobs)
        assert levels(f) == {"ok"}, f

    def test_running_job_blocks(self, case):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, lambda: {CASE})
        assert runstate.blocked(f)

    def test_no_squeue_blocks(self, case):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, lambda: None)
        assert runstate.blocked(f)

    def test_skip_slurm_is_a_warning(self, case):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, None)
        assert not runstate.blocked(f)
        assert "warn" in levels(f)

    def test_mixed_rpointers_block(self, case):
        run, archive = case
        (run / "rpointer.drv").write_text(f"{CASE}.cpl.r.0099-01-01-00000.nc\n")
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, no_jobs)
        assert any("disagree" in x.message for x in f if x.level == "block")

    def test_history_past_restart_blocks(self, case):
        run, archive = case
        (run / f"{CASE}.cam.h0.0101-01.nc").write_text("h0")
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, no_jobs)
        assert any("stale" in x.message for x in f if x.level == "block")

    def test_missing_archive_set_blocks(self, case, tmp_path):
        run, _ = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", tmp_path / "nope", no_jobs)
        assert any("rollback" in x.message for x in f if x.level == "block")

    def test_jumped_archive_copy_blocks(self, case):
        run, archive = case
        restart.apply_ice_jump(archive / "rest" / D0 / f"{CASE}.cice.r.{D0}.nc", 1.2)
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, no_jobs)
        assert any("not a pristine" in x.message for x in f if x.level == "block")

    @pytest.mark.parametrize("year,level", [(100, "ok"), (98, "warn"), (90, "block"),
                                            (101, "block")])
    def test_advice_staleness(self, case, year, level):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, no_jobs,
                               {"case": CASE, "model_year": year,
                                "schema_version": ADVICE_SCHEMA_VERSION})
        assert levels(f) >= {level}
        if level == "ok":
            assert not runstate.blocked(f)

    def test_advice_for_other_case_blocks(self, case):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, no_jobs,
                               {"case": "other", "model_year": 100,
                                "schema_version": ADVICE_SCHEMA_VERSION})
        assert runstate.blocked(f)

    def test_advice_missing_schema_version_blocks(self, case):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, no_jobs,
                               {"case": CASE, "model_year": 100})
        assert runstate.blocked(f)
        assert any("schema_version" in x.message for x in f if x.level == "block")

    def test_advice_unknown_schema_version_blocks(self, case):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, no_jobs,
                               {"case": CASE, "model_year": 100,
                                "schema_version": "999"})
        assert runstate.blocked(f)
        assert any("schema_version" in x.message for x in f if x.level == "block")

    def test_advice_missing_case_blocks(self, case):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", archive, no_jobs,
                               {"model_year": 100, "schema_version": ADVICE_SCHEMA_VERSION})
        assert runstate.blocked(f)

    def test_no_archive_without_override_blocks(self, case):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", None, no_jobs)
        assert any("no --archive given" in x.message for x in f if x.level == "block")

    def test_no_archive_with_override_and_full_local_set_passes(self, case):
        run, archive = case
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", None, no_jobs,
                               allow_no_archive=True)
        assert not runstate.blocked(f), f
        assert any("--allow-no-archive-rollback" in x.message for x in f
                  if x.level == "ok")

    def test_no_archive_with_override_but_missing_component_blocks(self, case):
        run, archive = case
        (run / f"{CASE}.cam.r.{D0}.nc").unlink()
        f = runstate.preflight(run / f"{CASE}.cice.r.{D0}.nc", None, no_jobs,
                               allow_no_archive=True)
        assert runstate.blocked(f)
        assert any("missing from" in x.message for x in f if x.level == "block")


def advance(run, archive):
    """Jump at D0, then simulate the post-jump segment and its st_archive."""
    cice = run / f"{CASE}.cice.r.{D0}.nc"
    restart.apply_ice_jump(cice, 1.5)
    write_set(run, D1, seed=5)
    for m in range(1, 13):
        (run / f"{CASE}.cam.h0.0101-{m:02d}.nc").write_text("h0")
    st_archive(run, archive, D1)
    # st_archive removes the older restart set from the run directory
    for p in list(run.iterdir()):
        if p.is_file() and D0 in p.name:
            p.unlink()
    return cice


class TestRollback:
    def test_from_archive(self, case):
        run, archive = case
        v0 = read(archive / "rest" / D0 / f"{CASE}.cice.r.{D0}.nc", "vicen")
        cice = advance(run, archive)
        plan = runstate.plan_rollback(run, archive, probe=no_jobs)
        assert plan.source == "archive"
        assert not runstate.blocked(plan.findings), plan.findings
        assert any("0101-" in d for d in plan.discarded)          # post-jump h0
        assert any(D1 in d for d in plan.discarded)                # post-jump restarts
        assert not any(D0 in d for d in plan.discarded)
        runstate.execute_rollback(plan, archive)
        for rp in ("rpointer.ice", "rpointer.drv", "rpointer.ocn", "rpointer.atm"):
            assert D0 in (run / rp).read_text() and D1 not in (run / rp).read_text()
        np.testing.assert_array_equal(read(cice, "vicen"), v0)
        assert restart.find_jump_logs(run) == []
        saved = run / restart.STATE_DIR / f"rollback-{D0}" / "rpointer.ice"
        assert D1 in saved.read_text()

    def test_from_rundir_when_not_archived(self, tmp_path):
        run = tmp_path / "run"
        run.mkdir()
        write_set(run, D0)
        v0 = read(run / f"{CASE}.cice.r.{D0}.nc", "vicen")
        restart.apply_ice_jump(run / f"{CASE}.cice.r.{D0}.nc", 1.5)
        write_set(run, D1, seed=5)            # DOUT_S off: old set stays in run/
        plan = runstate.plan_rollback(run, None, probe=no_jobs)
        assert plan.source == "rundir"
        assert not runstate.blocked(plan.findings), plan.findings
        runstate.execute_rollback(plan)
        assert D0 in (run / "rpointer.atm").read_text()
        assert D1 not in (run / "rpointer.atm").read_text()
        np.testing.assert_array_equal(read(run / f"{CASE}.cice.r.{D0}.nc", "vicen"), v0)

    def test_rundir_rollback_blocks_when_files_archived_away(self, case):
        run, archive = case
        advance(run, archive)
        plan = runstate.plan_rollback(run, None, probe=no_jobs)
        assert runstate.blocked(plan.findings)

    def test_running_job_blocks(self, case):
        run, archive = case
        advance(run, archive)
        plan = runstate.plan_rollback(run, archive, probe=lambda: {CASE})
        assert runstate.blocked(plan.findings)

    def test_no_active_jump(self, case):
        run, archive = case
        with pytest.raises(FileNotFoundError):
            runstate.plan_rollback(run, archive, probe=no_jobs)

    def test_cli(self, case):
        run, archive = case
        advance(run, archive)
        assert main(["rollback", "--rundir", str(run), "--archive", str(archive),
                     "--skip-slurm-check", "--dry-run"]) == 0
        assert D1 in (run / "rpointer.ice").read_text()           # dry run
        assert main(["rollback", "--rundir", str(run), "--archive", str(archive),
                     "--skip-slurm-check", "--yes"]) == 0
        assert D0 in (run / "rpointer.ice").read_text()


def test_cli_jump_preflight_end_to_end(case, tmp_path):
    run, archive = case
    advice = tmp_path / "advice.json"
    advice.write_text(json.dumps({"case": CASE, "ice_factor": 1.25, "model_year": 100,
                                  "schema_version": ADVICE_SCHEMA_VERSION}))
    assert main(["jump", "--rundir", str(run), "--advice", str(advice),
                 "--archive", str(archive), "--skip-slurm-check", "--yes"]) == 0
    log = restart.find_jump_logs(run)[0]
    assert json.loads(log.read_text())["jump_model_year"] == 101


def test_cli_jump_blocks_without_archive_or_override(case, tmp_path):
    run, archive = case
    advice = tmp_path / "advice.json"
    advice.write_text(json.dumps({"case": CASE, "ice_factor": 1.25, "model_year": 100,
                                  "schema_version": ADVICE_SCHEMA_VERSION}))
    assert main(["jump", "--rundir", str(run), "--advice", str(advice),
                 "--skip-slurm-check", "--yes"]) == 3
    assert restart.find_jump_logs(run) == []


def test_cli_jump_with_unknown_schema_version_blocks(case, tmp_path):
    run, archive = case
    advice = tmp_path / "advice.json"
    advice.write_text(json.dumps({"case": CASE, "ice_factor": 1.25, "model_year": 100,
                                  "schema_version": "999"}))
    assert main(["jump", "--rundir", str(run), "--advice", str(advice),
                 "--archive", str(archive), "--skip-slurm-check", "--yes"]) == 3
    assert restart.find_jump_logs(run) == []
