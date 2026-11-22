"""Hard per-test time limit, so no test can hang.

Every test gets 10 seconds by default (override with @pytest.mark.timeout(s)).
If it runs over, it is killed with SIGALRM and reported as a failure.
run_tests.sh also caps the whole suite with an outer `timeout`.
"""

import signal

import pytest

DEFAULT_TIMEOUT = 10


@pytest.fixture(autouse=True)
def _hard_timeout(request):
    marker = request.node.get_closest_marker("timeout")
    seconds = marker.args[0] if marker else DEFAULT_TIMEOUT

    def on_alarm(signum, frame):
        pytest.fail(f"test exceeded hard limit of {seconds}s", pytrace=False)

    old = signal.signal(signal.SIGALRM, on_alarm)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)
