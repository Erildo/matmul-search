"""Import schemes written by ternary_flip_graph (txt format) into our Scheme.

Format: first line "n m p rank", then three lines with rank*n*m, rank*m*p and
rank*n*p numbers (U, V, W blocks, row-major). The tool writes W for C
transposed (Brent-equation convention), so by default we transpose each W block.
Our verifier is the judge: `load_fg` raises if the result does not verify.
"""

from __future__ import annotations

from pathlib import Path

from .scheme import Scheme
from .verify import verify


def _blocks(nums, r, rows, cols):
    if len(nums) != r * rows * cols:
        raise ValueError(f"expected {r * rows * cols} numbers, got {len(nums)}")
    size = rows * cols
    return [
        [nums[k * size + i * cols: k * size + (i + 1) * cols] for i in range(rows)]
        for k in range(r)
    ]


def parse_fg(text: str, transpose_w: bool = True) -> Scheme:
    lines = [l.split() for l in text.strip().splitlines() if l.strip()]
    n, m, p, r = map(int, lines[0])
    u, v, w = ([int(x) for x in line] for line in lines[1:4])
    U = _blocks(u, r, n, m)
    V = _blocks(v, r, m, p)
    if transpose_w:
        Wt = _blocks(w, r, p, n)
        W = [[[blk[l][i] for l in range(p)] for i in range(n)] for blk in Wt]
    else:
        W = _blocks(w, r, n, p)
    return Scheme(n, m, p, U, V, W, {"source": "ternary_flip_graph"})


def load_fg(path: str | Path, mod: int | None = None) -> Scheme:
    """Parse and verify. Tries the transposed-W convention first."""
    text = Path(path).read_text()
    for tw in (True, False):
        s = parse_fg(text, transpose_w=tw)
        if verify(s, mod=mod).ok:
            s.meta["w_transposed"] = tw
            return s
    raise ValueError(f"{path}: does not verify in either W convention (mod={mod})")
