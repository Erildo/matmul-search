"""Our importer for ternary_flip_graph output must agree with our verifier."""

import pytest

from matmul.import_fg import load_fg, parse_fg
from matmul.verify import verify

# a real 2x2x2 rank-7 scheme written by ternary_flip_graph (seed 1)
FG_222 = """2 2 2 7
0 0 1 1 1 0 -1 0 0 1 1 0 0 1 0 0 0 1 0 -1 0 0 1 0 1 1 0 0
0 0 1 0 1 -1 0 0 0 1 -1 0 0 1 0 -1 0 0 1 -1 1 0 -1 0 0 1 0 0
0 1 0 1 1 0 0 0 -1 0 0 1 0 0 -1 -1 0 0 0 1 1 1 0 0 1 0 1 0
"""


def test_real_tool_output_verifies(tmp_path):
    f = tmp_path / "s.txt"
    f.write_text(FG_222)
    s = load_fg(f)
    assert (s.format, s.rank) == ("2x2x2", 7)
    assert s.meta["w_transposed"] is True


def test_corrupted_tool_output_is_rejected(tmp_path):
    lines = FG_222.splitlines()
    nums = lines[1].split()
    nums[2] = "-1" if nums[2] == "1" else "1"
    lines[1] = " ".join(nums)
    f = tmp_path / "bad.txt"
    f.write_text("\n".join(lines))
    with pytest.raises(ValueError, match="does not verify"):
        load_fg(f)


def test_wrong_number_count_rejected():
    with pytest.raises(ValueError, match="expected"):
        parse_fg("2 2 2 7\n0 1\n0 1\n0 1\n")


def test_untransposed_reading_is_wrong_for_this_tool():
    assert not verify(parse_fg(FG_222, transpose_w=False)).ok
