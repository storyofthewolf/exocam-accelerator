import numpy as np
import pytest

from exocam_accelerate.trends import (
    TrendSeries,
    build_area_weights,
    fit_curvature,
    fit_tendency,
    horizontal_mean,
    layer_mean_series,
    segment_slopes,
)


@pytest.fixture
def grid():
    lat = np.linspace(-87.5, 87.5, 36)
    lon = np.arange(0.0, 360.0, 5.0)
    return lon, lat


class TestAreaWeights:
    def test_weights_sum_to_one(self, grid):
        lon, lat = grid
        w = build_area_weights(lon, lat)
        assert w.shape == (len(lat), len(lon))
        assert np.isclose(w.sum(), 1.0)

    def test_equator_outweighs_pole(self, grid):
        lon, lat = grid
        w = build_area_weights(lon, lat)
        i_eq = np.argmin(np.abs(lat))
        assert w[i_eq, 0] > w[0, 0]

    def test_uniform_field_mean_is_value(self, grid):
        lon, lat = grid
        w = build_area_weights(lon, lat)
        field = np.full((len(lat), len(lon)), 3.25)
        assert np.isclose(horizontal_mean(field, w), 3.25)

    def test_analytic_sin2_lat_mean(self, grid):
        # area mean of sin^2(lat) over the sphere is 1/3
        lon, lat = grid
        w = build_area_weights(lon, lat)
        field = np.sin(np.deg2rad(lat))[:, None] ** 2 * np.ones(len(lon))
        assert np.isclose(horizontal_mean(field, w), 1.0 / 3.0, atol=1e-3)


class TestHorizontalMean:
    def test_per_layer_shapes(self, grid):
        lon, lat = grid
        w = build_area_weights(lon, lat)
        nlev, nt = 4, 7
        field3 = np.random.default_rng(0).normal(size=(nlev, len(lat), len(lon)))
        assert horizontal_mean(field3, w).shape == (nlev,)
        field4 = np.random.default_rng(1).normal(size=(nt, nlev, len(lat), len(lon)))
        assert horizontal_mean(field4, w).shape == (nt, nlev)

    def test_masked_cells_excluded(self, grid):
        lon, lat = grid
        w = build_area_weights(lon, lat)
        field = np.ma.masked_array(
            np.full((len(lat), len(lon)), 2.0),
            mask=np.zeros((len(lat), len(lon)), dtype=bool),
        )
        field[:, :10] = 99.0
        field.mask[:, :10] = True
        assert np.isclose(horizontal_mean(field, w), 2.0)

    def test_shape_mismatch_raises(self, grid):
        lon, lat = grid
        w = build_area_weights(lon, lat)
        with pytest.raises(ValueError):
            horizontal_mean(np.zeros((3, 4)), w)


class TestTrendSeries:
    def test_validation(self):
        with pytest.raises(ValueError):
            TrendSeries(times=[0.0], values=[1.0])
        with pytest.raises(ValueError):
            TrendSeries(times=[0.0, 0.0], values=[1.0, 2.0])  # not increasing
        with pytest.raises(ValueError):
            TrendSeries(times=[0.0, 1.0], values=[1.0, 2.0, 3.0])  # length

    def test_window_length_and_state_shape(self):
        s = TrendSeries(times=[0.0, 1.0, 2.0], values=np.zeros((3, 5)))
        assert s.window_length == 2.0
        assert s.state_shape == (5,)


class TestFits:
    def test_fit_tendency_recovers_slope(self):
        t = np.linspace(0, 10, 40)
        rng = np.random.default_rng(42)
        values = 250.0 - 0.8 * t + rng.normal(0, 0.01, size=t.shape)
        s = TrendSeries(times=t, values=values)
        assert np.isclose(fit_tendency(s), -0.8, atol=0.01)

    def test_fit_tendency_per_layer(self):
        t = np.linspace(0, 10, 30)
        slopes = np.array([-1.0, 0.0, 2.5])
        values = 100.0 + t[:, None] * slopes[None, :]
        s = TrendSeries(times=t, values=values)
        assert np.allclose(fit_tendency(s), slopes)

    def test_fit_curvature_on_parabola(self):
        t = np.linspace(0, 4, 25)
        values = 1.0 + 2.0 * t + 3.0 * t**2   # d2X/dt2 = 6
        s = TrendSeries(times=t, values=values)
        b, d2 = fit_curvature(s)
        assert np.isclose(d2, 6.0)
        # slope at window midpoint t=2 is 2 + 6*2 = 14
        assert np.isclose(b, 14.0)

    def test_curvature_needs_three_samples(self):
        s = TrendSeries(times=[0.0, 1.0], values=[1.0, 2.0])
        with pytest.raises(ValueError):
            fit_curvature(s)

    def test_segment_slopes_sign_flip(self):
        # V-shaped series: negative then positive slope
        t = np.linspace(0, 10, 20)
        values = np.abs(t - 5.0)
        s = TrendSeries(times=t, values=values)
        slopes = segment_slopes(s, 2)
        assert slopes[0] < 0 < slopes[1]

    def test_segment_slopes_too_few_samples(self):
        s = TrendSeries(times=[0.0, 1.0, 2.0], values=[0.0, 1.0, 2.0])
        with pytest.raises(ValueError):
            segment_slopes(s, 2)


class TestLayerMeanSeries:
    def test_builds_per_layer_series(self, grid):
        lon, lat = grid
        w = build_area_weights(lon, lat)
        nt, nlev = 6, 3
        t = np.arange(nt, dtype=float)
        # each layer drifts at its own rate, uniformly in space
        rates = np.array([0.1, -0.5, 0.0])
        fields = (
            270.0
            + t[:, None, None, None] * rates[None, :, None, None]
            + np.zeros((nt, nlev, len(lat), len(lon)))
        )
        s = layer_mean_series(t, fields, w)
        assert s.values.shape == (nt, nlev)
        assert np.allclose(fit_tendency(s), rates)
