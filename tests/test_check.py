import pytest
from synth import stefan_columns, truncate

from exocam_accelerate.advise import AdvisorConfig, advise
from exocam_accelerate.check import CheckConfig, Verdict, check_jump

JUMP_YEAR = 61          # restart 0061-01-01: data through year 60 advised on


def jump_log(factor=None):
    adv = advise(truncate(stefan_columns(), JUMP_YEAR - 1), "synth",
                 AdvisorConfig(max_ice_factor=1.5))
    assert adv.jump, adv.reasons
    return {"jump_model_year": JUMP_YEAR,
            "ice_factor": factor if factor is not None else adv.ice_factor,
            "advice": adv.to_dict()}


def run_after(years_after, **kw):
    log = jump_log()
    cols = stefan_columns(years=JUMP_YEAR - 1 + years_after, jump_year=JUMP_YEAR,
                          factor=log["ice_factor"], **kw)
    return check_jump(cols, log, start_year=1), log


def test_pass_when_run_follows_the_law():
    res, _ = run_after(8)
    assert res.verdict is Verdict.PASS, res.reasons
    assert abs(res.metrics["land_error"]) < 0.02
    assert abs(res.metrics["dN"]) < 0.05


def test_wait_before_any_post_jump_year():
    log = jump_log()
    res = check_jump(truncate(stefan_columns(), JUMP_YEAR - 1), log)
    assert res.verdict is Verdict.WAIT


def test_wait_during_settling():
    res, _ = run_after(2)
    assert res.verdict is Verdict.WAIT
    assert res.settled_years == 0


def test_wait_with_too_few_settled_years():
    res, _ = run_after(3)
    assert res.verdict is Verdict.WAIT
    assert res.settled_years == 1


def test_fail_when_jump_never_landed():
    log = jump_log()
    cols = stefan_columns(years=JUMP_YEAR + 5)          # no jump in the data
    res = check_jump(cols, log)
    assert res.verdict is Verdict.FAIL
    assert "did not land" in res.reasons[0]


def test_fail_when_off_the_conduction_law():
    res, _ = run_after(8, N_offset_after=1.0)
    assert res.verdict is Verdict.FAIL
    assert any("conduction law" in r for r in res.reasons)


def test_early_fail_on_large_deviation():
    res, _ = run_after(3, N_offset_after=3.0)
    assert res.verdict is Verdict.FAIL


def test_fail_when_ice_melts_back():
    res, _ = run_after(8, melt_after=True)
    assert res.verdict is Verdict.FAIL
    assert any("melting back" in r for r in res.reasons)


def test_exit_codes():
    assert [v.exit_code for v in Verdict] == [0, 10, 20]


def test_log_without_fit_is_rejected():
    with pytest.raises(ValueError, match="conduction-law fit"):
        check_jump(stefan_columns(), {"jump_model_year": 50, "ice_factor": 1.2,
                                      "advice": None})
