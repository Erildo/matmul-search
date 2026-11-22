"""A matrix multiplication scheme and its JSON file format.

Convention: C = A @ B with A n x m, B m x p, C n x p.
A scheme of rank r is r triples (U_k, V_k, W_k) with U_k n x m, V_k m x p,
W_k n x p. It computes

    M_k = (sum_ij U_k[i][j] * A[i][j]) * (sum_jl V_k[j][l] * B[j][l])
    C[i][l] = sum_k W_k[i][l] * M_k

Coefficients are exact rationals (ints, or strings like "1/2" in JSON).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path

Matrix = list[list[Fraction]]


def _to_fraction(x) -> Fraction:
    if isinstance(x, bool):
        raise ValueError(f"boolean is not a coefficient: {x!r}")
    if isinstance(x, float):
        raise ValueError(f"floats are not allowed (inexact): {x!r}; use an int or 'p/q' string")
    return Fraction(x)


def _matrix(rows, shape: tuple[int, int], what: str) -> Matrix:
    rows_n, cols_n = shape
    if len(rows) != rows_n or any(len(r) != cols_n for r in rows):
        raise ValueError(f"{what} must be {rows_n}x{cols_n}")
    return [[_to_fraction(x) for x in r] for r in rows]


@dataclass
class Scheme:
    n: int
    m: int
    p: int
    U: list[Matrix]
    V: list[Matrix]
    W: list[Matrix]
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        for name, val in (("n", self.n), ("m", self.m), ("p", self.p)):
            if not isinstance(val, int) or val < 1:
                raise ValueError(f"{name} must be a positive int, got {val!r}")
        if not (len(self.U) == len(self.V) == len(self.W)):
            raise ValueError(
                f"U, V, W must have the same length (got {len(self.U)}, {len(self.V)}, {len(self.W)})"
            )
        self.U = [_matrix(u, (self.n, self.m), f"U[{k}]") for k, u in enumerate(self.U)]
        self.V = [_matrix(v, (self.m, self.p), f"V[{k}]") for k, v in enumerate(self.V)]
        self.W = [_matrix(w, (self.n, self.p), f"W[{k}]") for k, w in enumerate(self.W)]

    @property
    def rank(self) -> int:
        return len(self.U)

    @property
    def format(self) -> str:
        return f"{self.n}x{self.m}x{self.p}"

    # ---- JSON -------------------------------------------------------------

    def to_dict(self) -> dict:
        def enc(mat):
            return [[x.numerator if x.denominator == 1 else str(x) for x in r] for r in mat]

        return {
            "n": self.n,
            "m": self.m,
            "p": self.p,
            "rank": self.rank,
            "meta": self.meta,
            "U": [enc(u) for u in self.U],
            "V": [enc(v) for v in self.V],
            "W": [enc(w) for w in self.W],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Scheme":
        s = cls(d["n"], d["m"], d["p"], d["U"], d["V"], d["W"], d.get("meta", {}))
        if "rank" in d and d["rank"] != s.rank:
            raise ValueError(f"declared rank {d['rank']} but file has {s.rank} products")
        return s

    def save(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(self.to_dict(), indent=1) + "\n")

    @classmethod
    def load(cls, path: str | Path) -> "Scheme":
        return cls.from_dict(json.loads(Path(path).read_text()))
