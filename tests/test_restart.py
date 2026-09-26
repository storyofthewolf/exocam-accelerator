import json

import numpy as np
import pytest

netCDF4 = pytest.importorskip("netCDF4")

from exocam_accelerate import restart  # noqa: E402
from exocam_accelerate.cli import main  # noqa: E402

NAME = "case.cice.r.0101-01-01-00000.nc"


def make_cice_r(path):
    rng = np.random.default_rng(0)
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
    make_cice_r(tmp_path / NAME)
    (tmp_path / "rpointer.ice").write_text(f"./{NAME}\n")
    return tmp_path


def test_locate_via_rpointer(rundir):
    assert restart.locate_cice_restart(rundir) == (rundir / NAME).resolve()


def test_backup_name():
    assert restart.backup_path("/x/" + NAME).name == \
        "case.cice.r.0101-01-01-00000.pre-accel.nc"


def test_jump_scales_ice_and_keeps_area(rundir):
    path = rundir / NAME
    v0, e0, a0, s0 = (read(path, n) for n in ("vicen", "eicen", "aicen", "vsnon"))
    rec = restart.apply_ice_jump(path, 1.5, provenance={"case": "case"})
    assert rec.written
    np.testing.assert_allclose(read(path, "vicen"), 1.5 * v0)
    np.testing.assert_allclose(read(path, "eicen"), 1.5 * e0)
    np.testing.assert_array_equal(read(path, "aicen"), a0)
    np.testing.assert_array_equal(read(path, "vsnon"), s0)
    # pristine backup holds the original
    np.testing.assert_array_equal(read(rec.backup, "vicen"), v0)
    # provenance recorded in file and sidecar
    with netCDF4.Dataset(path) as ds:
        meta = json.loads(ds.getncattr(restart.ATTR))
        assert ds.istep1 == 87600
    assert meta["ice_factor"] == 1.5
    log = json.loads(restart.log_path(path).read_text())
    assert log["advice"] == {"case": "case"}


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


def test_restore(rundir):
    path = rundir / NAME
    v0 = read(path, "vicen")
    restart.apply_ice_jump(path, 1.5)
    restart.restore(path)
    np.testing.assert_array_equal(read(path, "vicen"), v0)
    assert not restart.log_path(path).exists()


def test_snow_factor(rundir):
    path = rundir / NAME
    s0 = read(path, "vsnon")
    restart.apply_ice_jump(path, 1.5, snow_factor=1.1)
    np.testing.assert_allclose(read(path, "vsnon"), 1.1 * s0)


def test_cli_jump_from_advice(rundir, tmp_path, capsys):
    path = rundir / NAME
    v0 = read(path, "vicen")
    advice = tmp_path / "advice.json"
    advice.write_text(json.dumps({"case": "case", "ice_factor": 1.3}))
    assert main(["jump", "--rundir", str(rundir), "--advice", str(advice), "--yes"]) == 0
    np.testing.assert_allclose(read(path, "vicen"), 1.3 * v0)
    assert main(["restore", "--rundir", str(rundir)]) == 0
    np.testing.assert_array_equal(read(path, "vicen"), v0)


def test_cli_refuses_no_jump_advice(rundir, tmp_path):
    advice = tmp_path / "advice.json"
    advice.write_text(json.dumps({"case": "case", "ice_factor": None,
                                  "reasons": ["ice edge still moving"]}))
    assert main(["jump", "--rundir", str(rundir), "--advice", str(advice), "--yes"]) == 1


def test_cli_hard_bound(rundir):
    assert main(["jump", "--rundir", str(rundir), "--ice-factor", "3", "--yes"]) == 2
