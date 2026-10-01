import json
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

import numpy as np
import pytest
from synth import stefan_columns, truncate

from exocam_accelerate.advise import AdvisorConfig, advise
from exocam_accelerate.viewer import (case_payload, case_summary, config_from_query,
                                      find_jump_logs, make_handler, scan_directory)

JUMP_YEAR = 61


def write_case(directory, case, cols, years, start=1):
    """Write synthetic columns as exocam-trend cam/cice text files."""
    span = f"{start:04d}-01-{start + years - 1:04d}-12"
    cam = ["energy_top", "TS", "ICEFRAC"]
    cice = ["hi", "qi", "Tsfc"]
    for comp, names in (("cam", cam), ("cice", cice)):
        keys = [f"{n}_{w}" for n in names for w in ("native", "int2") if f"{n}_{w}" in cols]
        header = "month  " + "  ".join(keys)
        rows = np.column_stack([cols["month"]] + [cols[k] for k in keys])
        np.savetxt(directory / f"{case}_{span}_{comp}.txt", rows, header=header,
                   comments="", fmt="%.8g")


def jump_log(case="synth"):
    adv = advise(truncate(stefan_columns(), JUMP_YEAR - 1), case,
                 AdvisorConfig(max_ice_factor=1.5))
    return {"case": case, "jump_model_year": JUMP_YEAR, "restart_date":
            f"{JUMP_YEAR:04d}-01-01-00000", "ice_factor": adv.ice_factor,
            "advice": adv.to_dict()}


def test_payload_before_jump_is_json_safe_and_complete():
    cols = truncate(stefan_columns(), 60)
    p = case_payload(cols, "synth")
    json.dumps(p, allow_nan=False)
    assert p["advice"]["ice_factor"] is not None
    assert p["ice"]["after"]["hi"] > p["ice"]["now"]["hi"]
    assert len(p["ice"]["curve"]["hi"]) == len(p["ice"]["curve"]["N"])
    assert p["ice"]["stefan"]["slope_h2"] > 0
    assert p["window"][1] == 60
    assert p["temperatures"]["TS"]["accepted"]
    assert {"energy_top", "hi", "TS", "ICEFRAC", "qi_per_hi"} <= set(p["series"])
    assert p["jumps"] == []


def test_payload_scores_a_jump_with_check():
    log = jump_log()
    cols = stefan_columns(years=JUMP_YEAR + 7, jump_year=JUMP_YEAR,
                          factor=log["ice_factor"])
    p = case_payload(cols, "synth", jump_logs=[log])
    json.dumps(p, allow_nan=False)
    (j,) = p["jumps"]
    assert j["year"] == JUMP_YEAR
    assert j["check"]["verdict"] == "PASS"
    assert j["advice"]["hi_before"] > 0
    assert case_summary(p)["verdict"] == "PASS"


def test_refused_case_still_renders():
    cols = truncate(stefan_columns(icefrac_drift=0.5), 60)
    p = case_payload(cols, "synth")
    json.dumps(p, allow_nan=False)
    assert p["advice"]["ice_factor"] is None
    assert p["advice"]["reasons"]
    assert case_summary(p)["refused"]


def test_config_from_query():
    cfg = config_from_query({"n_fraction": ["0.3"], "max_ice_factor": ["1.2"],
                             "window_years": ["auto"], "which": ["int1"],
                             "since_year": [""]})
    assert cfg.n_fraction == 0.3 and cfg.max_ice_factor == 1.2
    assert cfg.window_years is None and cfg.which == "int1" and cfg.since_year is None
    with pytest.raises(ValueError):
        config_from_query({"max_ice_factor": ["3"]})


def test_scan_uses_latest_span_and_finds_logs(tmp_path):
    write_case(tmp_path, "caseA", truncate(stefan_columns(), 40), 40)
    write_case(tmp_path, "caseA", truncate(stefan_columns(), 60), 60)
    sub = tmp_path / "logs"
    sub.mkdir()
    (sub / "caseA.cice.r.0061-01-01-00000.nc.accel.json").write_text(
        json.dumps(jump_log("caseA")))
    (sub / "caseA.cice.r.0041-01-01-00000.nc.accel.json.restored").write_text("{}")
    info = scan_directory(tmp_path)
    assert list(info) == ["caseA"]
    assert info["caseA"]["end"] == "0060-12"
    assert len(info["caseA"]["files"]) == 2
    logs = find_jump_logs(tmp_path)
    assert [g["jump_model_year"] for g in logs["caseA"]] == [JUMP_YEAR]


@pytest.fixture
def server(tmp_path):
    write_case(tmp_path, "caseA", truncate(stefan_columns(), 60), 60)
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), make_handler(tmp_path))
    th = threading.Thread(target=httpd.serve_forever, daemon=True)
    th.start()
    yield f"http://127.0.0.1:{httpd.server_address[1]}"
    httpd.shutdown()
    httpd.server_close()


def get(url):
    try:
        with urllib.request.urlopen(url, timeout=10) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


def test_server_routes(server):
    code, body = get(server + "/")
    assert code == 200 and b"Spin-up Viewer" in body
    code, body = get(server + "/api/cases")
    rows = json.loads(body)["cases"]
    assert code == 200 and rows[0]["case"] == "caseA" and rows[0]["ice_factor"]
    code, body = get(server + "/api/case?case=caseA&n_fraction=0.3&max_ice_factor=1.2")
    p = json.loads(body)
    assert code == 200 and p["advice"]["ice_factor"] <= 1.2
    code, body = get(server + "/api/case?case=caseA&max_ice_factor=5")
    assert code == 400 and "max_ice_factor" in json.loads(body)["error"]
    code, _ = get(server + "/api/case?case=nope")
    assert code == 404
    code, body = get(server + "/api/stamp")
    assert code == 200 and json.loads(body)["stamp"] > 0


def test_ocean_regime_payload():
    from synth import gregory_columns, truncate
    from exocam_accelerate.viewer import any_payload, case_regime, case_summary
    cols = truncate(gregory_columns(), 30)
    assert case_regime(cols) == "ocean"
    p = any_payload(cols, "hot", 1, {"n_fraction": ["0.5"]}, None, [])
    assert p["regime"] == "ocean" and p["gregory"]["after"] is not None
    assert p["gregory"]["projection"]["TS_jumped"] is not None
    assert "gap" in p["series"] and p["advice"]["windows_tried"]
    assert p["gregory"]["local"] and p["gregory"]["alpha_diff"] > 0
    line = any_payload(cols, "hot", 1, {"n_fraction": ["0.5"]}, None, [])
    assert line["gregory"]["curve"]["N"] is not None
    row = case_summary(p)
    assert row["regime"] == "ocean" and row["somtp_dT"] > 0 and not row["refused"]
    import json
    json.dumps(p)                     # JSON-safe
