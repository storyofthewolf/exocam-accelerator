"""Synthetic exocam-trend columns obeying conduction-limited ice growth."""

import numpy as np

A, B, C0, C1 = 0.1, -70.0, 210.0, 1.5


def stefan_columns(years=100, a=A, b=B, h0=15.0, k=6.0, icefrac=0.8,
                   icefrac_drift=0.0, jump_year=None, factor=1.0,
                   N_offset_after=0.0, melt_after=False, noise=0.0, seed=0):
    """Monthly columns: hi(t) = sqrt(h0^2 + 2 k t), N = a + b/hi, TS linear in N.

    With ``jump_year`` (model year, series starting at year 1), hi is multiplied
    by ``factor`` from that year on and keeps growing by the same Stefan law;
    ``N_offset_after`` shifts N off the law after the jump (a failed jump);
    ``melt_after`` makes hi decline after the jump. int2 = native here.
    """
    rng = np.random.default_rng(seed)
    month = np.arange(1, 12 * years + 1, dtype=float)
    t = month / 12.0
    hi = np.sqrt(h0**2 + 2 * k * t)
    if jump_year is not None:
        after = t > (jump_year - 1)
        t_j = jump_year - 1
        h_j = np.sqrt(h0**2 + 2 * k * t_j) * factor
        hi = np.where(after, np.sqrt(h_j**2 + 2 * k * (t - t_j)), hi)
        if melt_after:
            hi = np.where(after, h_j - 0.5 * (t - t_j), hi)
    N = a + b / hi
    if jump_year is not None:
        N = np.where(t > (jump_year - 1), N + N_offset_after, N)
    N = N + noise * rng.standard_normal(N.size)
    series = {
        "hi": hi,
        "energy_top": N,
        "TS": C0 + C1 * N,
        "Tsfc": -60.0 + C1 * N,
        "qi": -3.0e20 * hi,
        "ICEFRAC": icefrac + icefrac_drift * t / years,
    }
    cols = {"month": month}
    for name, v in series.items():
        cols[f"{name}_native"] = v
        cols[f"{name}_int2"] = v
    return cols


def truncate(cols, years):
    return {k: v[: 12 * years] for k, v in cols.items()}
