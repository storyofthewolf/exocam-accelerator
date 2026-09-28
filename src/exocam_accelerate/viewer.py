"""Local interactive viewer: ``exocam-accelerate view DIR``.

Serves a single-page app on 127.0.0.1 that shows, for each case whose
exocam-trend files sit in DIR, what the advisor and the post-jump check see:

* the phase-space picture — N against hi with the fitted conduction law
  ``N = a + b/hi``, the current point, the target, and where a (clipped) jump
  lands; post-jump years plotted over it, so "did the run land on the law?"
  is visible at a glance;
* ``TS`` against N with the Gregory-style reference;
* time series of N, TS, hi, ICEFRAC and qi/hi with the fit window, jump
  markers, the predicted post-jump levels, and the unaccelerated Stefan
  projection of hi (h^2 linear in t) that a jump short-cuts;
* live "what-if" advice: the advisor re-runs on every control change
  (n_fraction, max ice factor, window, ...), so the effect of a setting is
  seen before anything is written.

Everything shown is computed by ``advise`` and ``check_jump`` — the viewer adds
no physics of its own, so what it draws is what ``advise``/``jump``/``check``
would decide.

Inputs, all read from DIR (and one level of subdirectories for jump logs):

* ``<case>_<YYYY>-<MM>-<YYYY>-<MM>_{cam,cice,clm}.txt`` — exocam-trend output;
  when several spans of one case exist, the one ending latest is used;
* ``*.accel.json`` — jump logs copied from ``run/exocam_accelerate/``; each
  becomes a jump marker and is scored with ``check_jump``.

The page re-fetches when any input file changes (a local mtime probe every
few seconds — no cluster contact; fetching new data from the HPC is a
separate step). Standard library only (``http.server``); the page loads
Plotly from a CDN.
"""

from __future__ import annotations

import json
import math
import re
import threading
import webbrowser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import parse_qs, urlparse

import numpy as np

from .advise import (ICE_AREA_VAR, ICE_ENTHALPY_VAR, ICE_VOLUME_VAR, AdvisorConfig,
                     _annual, _have, advise, model_years)
from .check import check_jump
from .phase_space import PhaseGateConfig
from .restart import LOG_SUFFIX
from .trend_io import merge_trend_files

_SPAN = re.compile(r"^(?P<case>.+)_(?P<y0>\d{4})-(?P<m0>\d{2})-(?P<y1>\d{4})-"
                   r"(?P<m1>\d{2})_(?P<comp>cam|cice|clm)\.txt$")

#: series drawn in the time-series panel (those present are sent)
SERIES = ("energy_top", "energy_bot", "TS", "Tsfc", ICE_VOLUME_VAR, "hs",
          ICE_AREA_VAR, ICE_ENTHALPY_VAR)


# --------------------------------------------------------------------------
# directory scan
# --------------------------------------------------------------------------

def scan_directory(directory) -> Dict[str, dict]:
    """``{case: {start_year, end, files}}`` for the latest-ending span per case."""
    spans: Dict[tuple, List[Path]] = {}
    for path in sorted(Path(directory).glob("*.txt")):
        m = _SPAN.match(path.name)
        if not m or m.group("m0") != "01":
            continue
        key = (m.group("case"), int(m.group("y0")),
               (int(m.group("y1")), int(m.group("m1"))))
        spans.setdefault(key, []).append(path)
    out: Dict[str, dict] = {}
    for (case, y0, end), files in sorted(spans.items(), key=lambda kv: kv[0][2]):
        out[case] = {"start_year": y0, "end": f"{end[0]:04d}-{end[1]:02d}",
                     "files": [str(f) for f in files]}
    return dict(sorted(out.items()))


def find_jump_logs(directory) -> Dict[str, List[dict]]:
    """Jump logs (``*.nc.accel.json``) in DIR and its immediate subdirectories."""
    d = Path(directory)
    found = list(d.glob(f"*.nc{LOG_SUFFIX}")) + list(d.glob(f"*/*.nc{LOG_SUFFIX}"))
    out: Dict[str, List[dict]] = {}
    for path in sorted(found):
        try:
            log = json.loads(path.read_text())
        except (OSError, ValueError):
            continue
        if "jump_model_year" in log and "case" in log:
            log["_file"] = path.name
            out.setdefault(log["case"], []).append(log)
    for logs in out.values():
        logs.sort(key=lambda g: int(g["jump_model_year"]))
    return out


def stamp(directory) -> float:
    """Newest mtime of any viewer input in DIR (drives the page's auto-refresh)."""
    d = Path(directory)
    paths = list(d.glob("*.txt")) + list(d.glob(f"*.nc{LOG_SUFFIX}")) + \
        list(d.glob(f"*/*.nc{LOG_SUFFIX}"))
    return max((p.stat().st_mtime for p in paths), default=0.0)


# --------------------------------------------------------------------------
# payload (pure: columns in, JSON-able dict out)
# --------------------------------------------------------------------------

def clean(x):
    """Recursively make numpy/NaN/inf JSON-safe (non-finite -> None)."""
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return [clean(v) for v in x.tolist()]
    if isinstance(x, (np.floating, float)):
        return float(x) if math.isfinite(float(x)) else None
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def _fit_params(ext, form):
    try:
        return tuple(ext.fits[form].params)
    except (AttributeError, KeyError):
        return None


def case_payload(columns: Dict[str, np.ndarray], case: str, start_year: int = 1,
                 config: AdvisorConfig = AdvisorConfig(),
                 jump_logs: Optional[List[dict]] = None,
                 provenance: Optional[Dict[str, str]] = None) -> dict:
    """Everything the page draws for one case, as a JSON-able dict."""
    t, _ = _annual(columns, config.imbalance, "native")
    years = model_years(t, start_year)

    series = {}
    for var in SERIES:
        if not _have(columns, var, "native"):
            continue
        entry = {"native": _annual(columns, var, "native")[1]}
        if f"{var}_{config.which}" in columns:
            entry["smooth"] = _annual(columns, var, config.which)[1]
        series[var] = entry
    if ICE_VOLUME_VAR in series and ICE_ENTHALPY_VAR in series:
        h, q = series[ICE_VOLUME_VAR], series[ICE_ENTHALPY_VAR]
        with np.errstate(divide="ignore", invalid="ignore"):
            series["qi_per_hi"] = {k: q[k] / h[k] for k in h if k in q}

    adv = advise(columns, case, config, start_year=start_year, provenance=provenance)
    advice = adv.to_dict()

    # the fit window, in model years (same rule as advise)
    w = t > t[-1] - adv.window_years
    if adv.since_year is not None:
        w &= years >= adv.since_year + config.settle_years
    window = [int(years[w][0]), int(years[w][-1])] if w.any() else None

    # conduction law, current point, target, landing point
    ice = None
    ab = _fit_params(adv.ice, "hyperbolic")
    if ab is not None and ICE_VOLUME_VAR in series:
        a, b = ab
        h_nat = series[ICE_VOLUME_VAR]["native"]
        h_now = float(h_nat[-1])
        h_target = adv.ice.prediction if adv.ice is not None else float("nan")
        h_after = h_now * adv.ice_factor if adv.ice_factor else float("nan")
        ends = [h_now] + [x for x in (h_target, h_after) if np.isfinite(x)]
        lo = float(np.nanmin(h_nat[w])) if w.any() else h_now
        hi_top = min(max(ends) * 1.15, h_now * 4.0)
        grid = np.linspace(max(lo * 0.85, 1e-3), hi_top, 200)
        ice = {"a": a, "b": b, "corr": adv.ice.fits["hyperbolic"].corr,
               "curve": {"hi": grid, "N": a + b / grid},
               "now": {"hi": h_now, "N": adv.N_now_fit},
               "target": {"hi": h_target, "N": adv.N_target},
               "after": ({"hi": h_after, "N": adv.N_after}
                         if adv.ice_factor else None)}
        # unaccelerated Stefan projection: h^2 linear in t over the fit window
        which = adv.which_used
        h_fit = _annual(columns, ICE_VOLUME_VAR, which)[1]
        if w.sum() >= 2:
            s = float(np.polyfit(t[w], h_fit[w] ** 2, 1)[0])
            if s > 0:
                span = max(20.0, 1.15 * (adv.years_skipped or 0.0))
                dt = np.linspace(0.0, span, 120)
                ice["stefan"] = {"slope_h2": s, "year": int(years[-1]) + 0.5 + dt,
                                 "hi": np.sqrt(h_now**2 + s * dt)}

    temps = {}
    for var, ext in adv.temperatures.items():
        c = _fit_params(ext, "linear")
        if c is None:
            continue
        Nw = np.asarray(series["energy_top"]["native"])[w] if "energy_top" in series \
            else np.array([])
        pts = [x for x in [adv.N_now_fit, adv.N_after, adv.N_target]
               if x is not None and np.isfinite(x)]
        span = np.concatenate([Nw, pts]) if len(pts) or Nw.size else np.array([0.0])
        ng = np.linspace(float(np.nanmin(span)) - 0.1, float(np.nanmax(span)) + 0.1, 50)
        temps[var] = {"c0": c[0], "c1": c[1], "accepted": ext.accepted,
                      "corr": ext.fits["linear"].corr,
                      "curve": {"N": ng, var: c[0] + c[1] * ng},
                      "prediction": ext.prediction}

    jumps = []
    for log in jump_logs or []:
        entry = {"file": log.get("_file"), "year": int(log["jump_model_year"]),
                 "restart_date": log.get("restart_date"),
                 "ice_factor": log.get("ice_factor"), "time_utc": log.get("time_utc")}
        la = log.get("advice") or {}
        entry["advice"] = {k: la.get(k) for k in
                           ("N_now_fit", "N_target", "N_after", "years_skipped",
                            "ice_factor_raw", "model_year", "expected_after_jump")}
        try:
            entry["advice"]["hi_before"] = la["ice"]["X_now"]
            entry["advice"]["law"] = la["ice"]["fits"]["hyperbolic"]["params"]
        except (KeyError, TypeError):
            pass
        try:
            res = check_jump(columns, log, start_year)
            entry["check"] = {"verdict": res.verdict.name, "reasons": list(res.reasons),
                              "post_years": res.years_after,
                              "settled_years": res.settled_years,
                              "metrics": dict(res.metrics)}
        except (KeyError, ValueError) as e:
            entry["check"] = {"verdict": None, "reasons": [str(e)]}
        jumps.append(entry)

    return clean({"case": case, "start_year": start_year, "years": years,
                  "which": adv.which_used, "series": series, "advice": advice,
                  "window": window, "ice": ice, "temperatures": temps,
                  "jumps": jumps})


def case_summary(payload: dict) -> dict:
    """One row of the case list: advice headline and latest check verdict."""
    adv = payload["advice"]
    last = payload["jumps"][-1] if payload["jumps"] else None
    return {"case": payload["case"], "model_year": adv["model_year"],
            "N_now": adv["N_now"], "ice_factor": adv["ice_factor"],
            "refused": adv["ice_factor"] is None, "n_reasons": len(adv["reasons"]),
            "since_year": adv["since_year"], "n_jumps": len(payload["jumps"]),
            "verdict": (last or {}).get("check", {}).get("verdict") if last else None}


# --------------------------------------------------------------------------
# request -> config
# --------------------------------------------------------------------------

def _num(q, key, cast=float):
    v = q.get(key, [""])[0].strip()
    return cast(v) if v not in ("", "auto", "none", "null") else None


def config_from_query(q: Dict[str, List[str]]) -> AdvisorConfig:
    kw = {}
    for key, cast in (("n_fraction", float), ("max_ice_factor", float),
                      ("window_years", float), ("N_target", float),
                      ("since_year", int)):
        v = _num(q, key, cast)
        if v is not None:
            kw[key] = v
    which = q.get("which", [""])[0]
    if which in ("native", "int1", "int2"):
        kw["which"] = which
    ratio = _num(q, "max_extrapolation_ratio")
    if ratio is not None:
        kw["gate"] = PhaseGateConfig(max_extrapolation_ratio=ratio)
    return AdvisorConfig(**kw)


# --------------------------------------------------------------------------
# server
# --------------------------------------------------------------------------

def _page() -> bytes:
    return (resources.files("exocam_accelerate") / "viewer_assets" / "index.html"
            ).read_bytes()


def make_handler(directory):
    directory = Path(directory)

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, fmt, *args):      # keep the terminal quiet
            pass

        def _send(self, code, body: bytes, ctype="application/json"):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def _json(self, obj, code=200):
            self._send(code, json.dumps(obj).encode())

        def do_GET(self):
            url = urlparse(self.path)
            q = parse_qs(url.query)
            try:
                if url.path in ("/", "/index.html"):
                    return self._send(200, _page(), "text/html; charset=utf-8")
                if url.path == "/api/stamp":
                    return self._json({"stamp": stamp(directory)})
                if url.path == "/api/cases":
                    cfg = config_from_query(q)
                    logs = find_jump_logs(directory)
                    rows = []
                    for case, info in scan_directory(directory).items():
                        try:
                            cols = merge_trend_files(info["files"], case)
                            p = case_payload(cols, case, info["start_year"], cfg,
                                             logs.get(case))
                            rows.append(dict(case_summary(p), end=info["end"]))
                        except Exception as e:           # one bad case != no list
                            rows.append({"case": case, "error": str(e)})
                    return self._json({"directory": str(directory.resolve()),
                                       "cases": rows})
                if url.path == "/api/case":
                    case = q.get("case", [""])[0]
                    info = scan_directory(directory).get(case)
                    if info is None:
                        return self._json({"error": f"no trend files for {case!r}"}, 404)
                    cfg = config_from_query(q)
                    cols = merge_trend_files(info["files"], case)
                    p = case_payload(cols, case, info["start_year"], cfg,
                                     find_jump_logs(directory).get(case))
                    p["files"] = [Path(f).name for f in info["files"]]
                    return self._json(p)
                return self._json({"error": "not found"}, 404)
            except ValueError as e:          # e.g. an out-of-domain control value
                return self._json({"error": str(e)}, 400)
            except Exception as e:           # pragma: no cover - surfaced in the page
                return self._json({"error": f"{type(e).__name__}: {e}"}, 500)

    return Handler


def serve(directory, host: str = "127.0.0.1", port: int = 8765,
          open_browser: bool = True) -> ThreadingHTTPServer:
    """Start the viewer (blocking). Binds to localhost only by default."""
    directory = Path(directory)
    if not directory.is_dir():
        raise FileNotFoundError(f"{directory} is not a directory")
    httpd = None
    for p in range(port, port + 20):          # busy port -> try the next ones
        try:
            httpd = ThreadingHTTPServer((host, p), make_handler(directory))
            break
        except OSError:
            continue
    if httpd is None:
        raise RuntimeError(f"no free port in {port}..{port + 19} on {host}")
    url = f"http://{host}:{httpd.server_address[1]}/"
    n = len(scan_directory(directory))
    print(f"exocam-accelerate viewer: {n} case(s) in {directory.resolve()}")
    print(f"  {url}   (Ctrl-C to stop)")
    if open_browser:
        threading.Timer(0.5, webbrowser.open, args=(url,)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()
    return httpd
