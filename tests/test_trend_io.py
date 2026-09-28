import hashlib

import numpy as np
import pytest

from exocam_accelerate.trend_io import file_provenance, load_case, read_trend_text, trend_series

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


# exocam-trend writes one file per component; a case merges them.
CAM_FIXTURE = """month  TS_native  TS_int1  TS_int2
1  250.0  250.0  250.0
2  250.5  250.2  250.2
"""
CICE_FIXTURE = """month  hi_native  hi_int1  hi_int2  vicen005_native  vicen005_int1  vicen005_int2
1  1.20  1.20  1.20  0.30  0.30  0.30
2  1.25  1.22  1.22  0.31  0.30  0.30
"""


@pytest.fixture
def case_dir(tmp_path):
    (tmp_path / "coldA_0001-01-0150-12_cam.txt").write_text(CAM_FIXTURE)
    (tmp_path / "coldA_0001-01-0150-12_cice.txt").write_text(CICE_FIXTURE)
    return tmp_path


class TestLoadCase:
    def test_merges_components_under_one_month_axis(self, case_dir):
        cols = load_case(case_dir, "coldA")
        # atmosphere TS and sea-ice hi/vicen live together after merge
        assert "TS_native" in cols and "hi_native" in cols and "vicen005_native" in cols
        assert np.array_equal(cols["month"], [1, 2])

    def test_variable_resolves_regardless_of_component(self, case_dir):
        cols = load_case(case_dir, "coldA")
        s = trend_series(cols, "hi", which="native")   # from the cice file
        assert np.isclose(s.values[-1], 1.25)

    def test_missing_case_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            load_case(tmp_path, "nope")

    def test_mismatched_month_axis_raises(self, tmp_path):
        (tmp_path / "b_0001-01-0002-12_cam.txt").write_text(CAM_FIXTURE)
        (tmp_path / "b_0001-01-0002-12_cice.txt").write_text(
            "month  hi_native  hi_int1  hi_int2\n1  1.0  1.0  1.0\n"
        )  # only one month -> axis differs
        with pytest.raises(ValueError):
            load_case(tmp_path, "b")


class TestFileProvenance:
    def test_hashes_every_component_file(self, case_dir):
        prov = file_provenance(case_dir, "coldA")
        assert set(prov) == {"coldA_0001-01-0150-12_cam.txt",
                             "coldA_0001-01-0150-12_cice.txt"}
        assert prov["coldA_0001-01-0150-12_cam.txt"] == hashlib.sha256(
            CAM_FIXTURE.encode()).hexdigest()

    def test_no_files_is_empty(self, tmp_path):
        assert file_provenance(tmp_path, "nope") == {}
