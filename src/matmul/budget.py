"""Hard limits for every search we run.

Rule of this repo: no search loop runs without a Budget. A budget caps wall
time and/or step count, and the search must call `tick()` every iteration.
When the budget runs out, `BudgetExceeded` is raised, and the search reports
what it has so far and exits. Nothing runs "until it finds something".

    budget = Budget(max_seconds=60, max_steps=1_000_000)
    try:
        while True:
            budget.tick()
            ...
    except BudgetExceeded:
        pass
    print(budget.report())
"""

from __future__ import annotations

import time


class BudgetExceeded(Exception):
    pass


class Budget:
    def __init__(self, max_seconds: float | None = None, max_steps: int | None = None):
        if max_seconds is None and max_steps is None:
            raise ValueError("a Budget needs max_seconds or max_steps (unbounded runs are not allowed)")
        if max_seconds is not None and max_seconds <= 0:
            raise ValueError("max_seconds must be > 0")
        if max_steps is not None and max_steps <= 0:
            raise ValueError("max_steps must be > 0")
        self.max_seconds = max_seconds
        self.max_steps = max_steps
        self.steps = 0
        self.start = time.monotonic()
        self.stop_reason: str | None = None

    @property
    def elapsed(self) -> float:
        return time.monotonic() - self.start

    def tick(self, n: int = 1) -> None:
        self.steps += n
        if self.max_steps is not None and self.steps > self.max_steps:
            self.stop_reason = "max_steps"
            raise BudgetExceeded(f"step budget of {self.max_steps} used up")
        if self.max_seconds is not None and self.elapsed > self.max_seconds:
            self.stop_reason = "max_seconds"
            raise BudgetExceeded(f"time budget of {self.max_seconds}s used up")

    def report(self) -> dict:
        return {
            "steps": self.steps,
            "seconds": round(self.elapsed, 3),
            "max_steps": self.max_steps,
            "max_seconds": self.max_seconds,
            "stop_reason": self.stop_reason or "finished",
        }
