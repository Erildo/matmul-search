"""The verifier must accept every correct scheme and reject every broken one.

If the verifier is wrong, every 'result' we produce is worthless, so it is
tested from both sides: known-good schemes pass, deliberately broken ones fail.
"""

import copy
import json
from fractions import Fraction
from pathlib import Path

import pytest

from matmul.known import kron, naive, strassen
from matmul.scheme import Scheme
from matmul.verify import apply_scheme, main, random_check, verify

SCHEMES_DIR = Path(__file__).resolve().parent.parent / "schemes"


# ---- known-good schemes pass ---------------------------------------------

@pytest.mark.parametrize("fmt", [(1, 1, 1), (2, 2, 2), (2, 3, 4), (3, 3, 3), (4, 2, 3)])
def test_naive_is_correct(fmt):
    s = naive(*fmt)
    assert s.rank == fmt[0] * fmt[1] * fmt[2]
    r = verify(s)
    assert r.ok, r.examples
    assert random_check(s)


def test_strassen_is_correct_over_Q():
    s = strassen()
    assert s.rank == 7
    assert verify(s).ok
    assert random_check(s, trials=5)


@pytest.mark.parametrize("prime", [2, 3, 5, 7])
def test_strassen_is_correct_mod_p(prime):
    assert verify(strassen(), mod=prime).ok


def test_strassen_squared_is_4x4_rank_49():
    s = kron(strassen(), strassen())
    assert (s.n, s.m, s.p, s.rank) == (4, 4, 4, 49)
    assert verify(s).ok
    assert random_check(s, trials=2)


def test_kron_of_rectangular_schemes():
    s = kron(naive(1, 2, 3), strassen())
    assert (s.format, s.rank) == ("2x4x6", 42)
    assert verify(s).ok


# ---- broken schemes fail ---------------------------------------------------

def _broken(mutate):
    s = strassen()
    d = copy.deepcopy(s.to_dict())
    mutate(d)
    return Scheme.from_dict(d)


def test_dropping_a_product_fails():
    s = strassen()
    bad = Scheme(2, 2, 2, s.U[:-1], s.V[:-1], s.W[:-1])
    r = verify(bad)
    assert not r.ok and r.wrong_entries > 0
    assert not random_check(bad)


@pytest.mark.parametrize("which", ["U", "V", "W"])
@pytest.mark.parametrize("k", range(7))
def test_flipping_any_single_sign_fails(which, k):
    # flip the first nonzero coefficient of product k in U, V or W
    def mutate(d):
        mat = d[which][k]
        for row in mat:
            for c, x in enumerate(row):
                if x != 0:
                    row[c] = -x
                    return
    bad = _broken(mutate)
    assert not verify(bad).ok
    assert not verify(bad, mod=3).ok


def test_correct_mod_2_but_not_over_Q():
    # multiplying a coefficient by 3 is invisible mod 2 but breaks it over Q
    def mutate(d):
        d["W"][0][0][0] = 3
    bad = _broken(mutate)
    assert verify(bad, mod=2).ok
    assert not verify(bad).ok
    assert not verify(bad, mod=3).ok


def test_swapping_products_between_slots_fails():
    def mutate(d):
        d["W"][0], d["W"][1] = d["W"][1], d["W"][0]
    assert not verify(_broken(mutate)).ok


def test_failure_report_names_entries():
    s = strassen()
    bad = Scheme(2, 2, 2, s.U[1:], s.V[1:], s.W[1:])
    r = verify(bad, max_examples=3)
    assert len(r.examples) == 3
    assert "WRONG" in r.summary()


# ---- rational coefficients -------------------------------------------------

def test_rational_coefficients_are_exact():
    # split a11*b11 into two halves: (1/2 a11) b11 + (1/2 a11) b11
    U = [[[Fraction(1, 2)]], [["1/2"]]]
    V = [[[1]], [[1]]]
    W = [[[1]], [[1]]]
    s = Scheme(1, 1, 1, U, V, W)
    assert verify(s).ok
    assert verify(s, mod=3).ok


def test_denominator_not_invertible_mod_p_raises():
    s = Scheme(1, 1, 1, [[["1/2"]], [["1/2"]]], [[[1]], [[1]]], [[[1]], [[1]]])
    with pytest.raises(ValueError, match="denominator"):
        verify(s, mod=2)


# ---- input validation ------------------------------------------------------

def test_floats_rejected():
    with pytest.raises(ValueError, match="float"):
        Scheme(1, 1, 1, [[[0.5]]], [[[1]]], [[[1]]])


def test_wrong_shapes_rejected():
    with pytest.raises(ValueError, match="U\\[0\\]"):
        Scheme(2, 2, 2, [[[1, 0]]], [[[1, 0], [0, 0]]], [[[1, 0], [0, 0]]])


def test_mismatched_counts_rejected():
    with pytest.raises(ValueError, match="same length"):
        Scheme(1, 1, 1, [[[1]]], [], [[[1]]])


def test_declared_rank_must_match(tmp_path):
    d = strassen().to_dict()
    d["rank"] = 6
    with pytest.raises(ValueError, match="declared rank"):
        Scheme.from_dict(d)


# ---- file round trip and CLI ----------------------------------------------

def test_save_load_roundtrip(tmp_path):
    s = kron(strassen(), naive(1, 1, 2))
    f = tmp_path / "s.json"
    s.save(f)
    t = Scheme.load(f)
    assert t.to_dict() == s.to_dict()
    assert verify(t).ok


def test_cli_exit_codes(tmp_path, capsys):
    good = tmp_path / "good.json"
    strassen().save(good)
    d = strassen().to_dict()
    d["U"][0][0][0] = 0
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps(d))
    assert main([str(good)]) == 0
    assert main([str(bad)]) == 1
    assert main([str(good), "--mod", "2"]) == 0


def test_apply_scheme_matches_hand_computation():
    A = [[1, 2], [3, 4]]
    B = [[5, 6], [7, 8]]
    assert apply_scheme(strassen(), A, B) == [[19, 22], [43, 50]]


# ---- every scheme committed to the repo must verify ----------------------

@pytest.mark.timeout(60)
@pytest.mark.parametrize("path", sorted(SCHEMES_DIR.glob("*.json")), ids=lambda p: p.name)
def test_every_committed_scheme_verifies(path):
    s = Scheme.load(path)
    field = s.meta.get("field", "Q")
    mod = None if field == "Q" else int(field.removeprefix("GF(").removesuffix(")"))
    r = verify(s, mod=mod)
    assert r.ok, r.summary()
    if "rank" in path.name:
        assert f"rank{s.rank}" in path.name, "file name must state the true rank"
