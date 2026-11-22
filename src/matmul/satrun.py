"""Run one cyclic F3 SAT instance: solve, decode, verify, save, and for UNSAT re-solve with a DRAT proof
and check it with drat-trim. Used by scripts/run_006.py and scripts/run_007.py.
"""

from __future__ import annotations

import subprocess
import time
from pathlib import Path

from .cyclic import expand, is_cyclic_invariant
from .verify import verify


def solver_cmd(solver, seed, cnf, proof=None):
    cmd = [solver, "-q", f"--seed={seed}", str(cnf)]
    if proof:
        cmd.append(str(proof))
    return cmd


def solve_instance(enc, name, solver, seed, cap, *, work: Path, found: Path, root: Path, drat_trim: Path,
                   attempt: str, prefix: str, check_cap=None, meta=None) -> dict:
    """Solve the CNF of `enc` (a CyclicF3CNF). Returns the result record (status, seconds, sizes, proof verdict)."""
    n = enc.n
    work.mkdir(parents=True, exist_ok=True)
    cnf = work / f"{name}.cnf"
    cnf.write_text(enc.dimacs())
    out = dict(meta or {})
    out.update(n=n, a=enc.a, b=enc.b, rank_max=3 * enc.a + enc.b, solver=solver, seed=seed, cap_seconds=cap,
               vars=enc.nv, clauses=len(enc.clauses), equations=len(enc.equations))
    t0 = time.time()
    try:
        r = subprocess.run(solver_cmd(solver, seed, cnf), capture_output=True, text=True, timeout=cap)
        status = {10: "SAT", 20: "UNSAT"}.get(r.returncode, f"exit{r.returncode}")
    except subprocess.TimeoutExpired:
        status = "TIMEOUT"
    out["seconds"] = round(time.time() - t0, 2)
    out["status"] = status
    if status == "SAT":
        model = [int(x) for line in r.stdout.splitlines() if line.startswith("v") for x in line.split()[1:]]
        orbits, fixed = enc.decode(model)
        s = expand(n, orbits, fixed, {"attempt": attempt, "instance": name, "field": "GF(3)"})
        out.update(found_rank=s.rank, verified_mod3=verify(s, mod=3).ok, cyclic_invariant=is_cyclic_invariant(s),
                   orbits_used=sum(1 for o in orbits if any(e for m in o for r_ in m for e in r_)),
                   fixed_used=sum(1 for x in fixed if any(e for r_ in x for e in r_)))
        found.mkdir(parents=True, exist_ok=True)
        path = found / f"{prefix}_{n}x{n}x{n}_rank{s.rank}_{name}.json"
        s.save(path)
        out["scheme_file"] = str(path.relative_to(root))
    if status == "UNSAT" and check_cap:
        proof = work / f"{name}.drat"
        t1 = time.time()
        try:
            r2 = subprocess.run(solver_cmd(solver, seed, cnf, proof), capture_output=True, text=True, timeout=cap * 2)
            out["proof_solve_status"] = {10: "SAT", 20: "UNSAT"}.get(r2.returncode, f"exit{r2.returncode}")
        except subprocess.TimeoutExpired:
            out["proof_solve_status"] = "TIMEOUT"
        out["proof_bytes"] = proof.stat().st_size if proof.exists() else 0
        if out["proof_solve_status"] == "UNSAT":
            try:
                r3 = subprocess.run([str(drat_trim), str(cnf), str(proof), "-t", str(int(check_cap))],
                                    capture_output=True, text=True, timeout=check_cap + 30)
                out["drat_trim"] = "VERIFIED" if "s VERIFIED" in r3.stdout else (
                    [l for l in r3.stdout.splitlines() if l.startswith("s ")] or ["no verdict"])[0]
            except subprocess.TimeoutExpired:
                out["drat_trim"] = "TIMEOUT"
        out["proof_seconds"] = round(time.time() - t1, 2)
        if out.get("drat_trim") == "VERIFIED":
            proof.unlink()  # large; reproducible from the CNF and seed
    return out
