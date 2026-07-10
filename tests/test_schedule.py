import pytest

from exocam_accelerate.schedule import DtSchedule


class TestValidation:
    def test_empty_raises(self):
        with pytest.raises(ValueError):
            DtSchedule(())

    def test_increasing_dt_raises(self):
        with pytest.raises(ValueError):
            DtSchedule(((5, 10.0), (5, 100.0)))

    def test_negative_dt_raises(self):
        with pytest.raises(ValueError):
            DtSchedule(((5, -1.0),))

    def test_zero_dt_only_final_phase(self):
        with pytest.raises(ValueError):
            DtSchedule(((5, 0.0), (5, 0.0)))
        DtSchedule(((5, 10.0), (5, 0.0)))  # trailing zero is fine (Turbet)

    def test_zero_iterations_raises(self):
        with pytest.raises(ValueError):
            DtSchedule(((0, 10.0),))


class TestLookup:
    def test_wordsworth2013(self):
        s = DtSchedule.wordsworth2013()
        assert s.total_iterations == 20
        assert s.dt_for(0) == 100.0
        assert s.dt_for(4) == 100.0
        assert s.dt_for(5) == 10.0
        assert s.dt_for(19) == 10.0

    def test_past_end_raises(self):
        s = DtSchedule.wordsworth2013()
        with pytest.raises(IndexError):
            s.dt_for(20)
        with pytest.raises(IndexError):
            s.dt_for(-1)

    def test_iteration_matches_lookup(self):
        s = DtSchedule.from_phases([(2, 50.0), (3, 10.0), (1, 0.0)])
        assert list(s) == [50.0, 50.0, 10.0, 10.0, 10.0, 0.0]
        assert [s.dt_for(i) for i in range(6)] == list(s)
