"""Budgets must actually stop runaway loops; run records must be complete."""

import time

import pytest

from matmul.budget import Budget, BudgetExceeded
from matmul.record import read_runs, record_run


def test_unbounded_budget_is_refused():
    with pytest.raises(ValueError, match="unbounded"):
        Budget()


@pytest.mark.parametrize("kw", [{"max_seconds": 0}, {"max_steps": 0}, {"max_seconds": -1}])
def test_nonpositive_limits_refused(kw):
    with pytest.raises(ValueError):
        Budget(**kw)


def test_step_budget_stops_infinite_loop():
    b = Budget(max_steps=1000)
    with pytest.raises(BudgetExceeded):
        while True:
            b.tick()
    assert b.steps == 1001
    assert b.report()["stop_reason"] == "max_steps"


@pytest.mark.timeout(5)
def test_time_budget_stops_infinite_loop():
    b = Budget(max_seconds=0.2)
    t0 = time.monotonic()
    with pytest.raises(BudgetExceeded):
        while True:
            b.tick()
    assert time.monotonic() - t0 < 1.0
    assert b.report()["stop_reason"] == "max_seconds"


def test_finished_run_reports_finished():
    b = Budget(max_steps=10)
    for _ in range(5):
        b.tick()
    assert b.report()["stop_reason"] == "finished"


def test_record_run_appends_and_reads(tmp_path):
    f = tmp_path / "runs.jsonl"
    base = dict(format="2x2x2", method="test", field="Q", best_rank_found=7,
                verified=True, budget={"max_steps": 10})
    record_run(f, attempt="000-a", **base)
    record_run(f, attempt="000-b", **base)
    runs = read_runs(f)
    assert [r["attempt"] for r in runs] == ["000-a", "000-b"]
    assert all("time" in r for r in runs)


def test_record_run_requires_fields(tmp_path):
    with pytest.raises(ValueError, match="missing"):
        record_run(tmp_path / "runs.jsonl", attempt="x")
