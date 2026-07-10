import numpy as np
import pytest

from exocam_accelerate.trend_io import read_trend_text, trend_series

# Matches the documented exocam-trend data/*.txt format:
# month  VAR_native  VAR_int1  VAR_int2 ... per variable
FIXTURE = """month  TS_native  TS_int1  TS_int2  ICEFRAC_native  ICEFRAC_int1  ICEFRAC_int2
1  284.28  284.41  284.41  0.0470  0.0474  0.0474
2  285.08  284.63  284.63  0.0487  0.0478  0.0478
3  286.65  285.14  285.14  0.0521  0.0489  0.0489
4  288.00  285.71  285.71  0.0542  0.0500  0.0500
"""


@pytest.fixture
def trend_file(tmp_path):
    path = tmp_path / "case_0001-01-0001-04_cam.txt"
    path.write_text(FIXTURE)
    return path


class TestReadTrendText:
    def test_parses_all_columns(self, trend_file):
        cols = read_trend_text(trend_file)
        assert set(cols) == {
            "month",
            "TS_native", "TS_int1", "TS_int2",
            "ICEFRAC_native", "ICEFRAC_int1", "ICEFRAC_int2",
        }
        assert np.array_equal(cols["month"], [1, 2, 3, 4])
        assert np.isclose(cols["TS_native"][2], 286.65)

    def test_rejects_non_trend_file(self, tmp_path):
        bad = tmp_path / "junk.txt"
        bad.write_text("time  TS\n1  284.0\n")
        with pytest.raises(ValueError):
            read_trend_text(bad)


class TestTrendSeries:
    def test_builds_series_in_years(self, trend_file):
        cols = read_trend_text(trend_file)
        s = trend_series(cols, "TS")
        assert np.allclose(s.times, np.array([1, 2, 3, 4]) / 12.0)
        assert s.values.shape == (4,)

    def test_which_and_window(self, trend_file):
        cols = read_trend_text(trend_file)
        s = trend_series(cols, "ICEFRAC", which="int1", last_n_months=2)
        assert len(s.times) == 2
        assert np.isclose(s.values[-1], 0.0500)

    def test_missing_variable_raises(self, trend_file):
        cols = read_trend_text(trend_file)
        with pytest.raises(KeyError):
            trend_series(cols, "FLNT")
