"""Attempt 005: SAT search for cyclic-invariant schemes over F3 (see attempts/005-cyclic-f3-sat.md).

Each instance (n, a, b) is encoded by matmul.cyclic.CyclicF3CNF, written to scratch/005/, and solved by
cadical or kissat (1 thread) under a hard timeout. SAT models are decoded, expanded, verified mod 3 and
checked for cyclic invariance; the scheme is saved to results/f3/. UNSAT instances are re-run with a DRAT
proof, which drat-trim checks. Every instance gets one line in results/runs.jsonl.

usage: python3 -I scripts/run_005.py STAGE        (STAGE = 0, A or B)
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.cyclic import CyclicF3CNF, expand, is_cyclic_invariant  # noqa: E402
from matmul.record import record_run  # noqa: E402
from matmul.verify import verify  # noqa: E402

WORK = ROOT / "scratch" / "005"
FOUND = ROOT / "results" / "f3"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
ATTEMPT = "005-cyclic-f3-sat"


def solver_cmd(solver, seed, cnf, proof=None):
    if solver == "cadical":
        cmd = ["cadical", "-q", f"--seed={seed}", str(cnf)]
    else:
        cmd = ["kissat", "-q", f"--seed={seed}", str(cnf)]
    if proof:
        cmd.append(str(proof))
    return cmd


def run_instance(n, a, b, solver, seed, cap, check_cap=None):
    name = f"n{n}_a{a}_b{b}_{solver}_s{seed}"
    WORK.mkdir(parents=True, exist_ok=True)
    enc = CyclicF3CNF(n, a, b)
    cnf = WORK / f"{name}.cnf"
    cnf.write_text(enc.dimacs())
    out = dict(n=n, a=a, b=b, rank_max=3 * a + b, solver=solver, seed=seed, cap_seconds=cap,
               vars=enc.nv, clauses=len(enc.clauses), equations=len(enc.equations))
    t0 = time.time()
    try:
        r = subprocess.run(solver_cmd(solver, seed, cnf), capture_output=True, text=True, timeout=cap)
        out["seconds"] = round(time.time() - t0, 2)
        status = {10: "SAT", 20: "UNSAT"}.get(r.returncode, f"exit{r.returncode}")
    except subprocess.TimeoutExpired:
        out["seconds"] = round(time.time() - t0, 2)
        status = "TIMEOUT"
    out["status"] = status
    if status == "SAT":
        model = [int(x) for line in r.stdout.splitlines() if line.startswith("v") for x in line.split()[1:]]
        orbits, fixed = enc.decode(model)
        s = expand(n, orbits, fixed, {"attempt": ATTEMPT, "instance": name, "field": "GF(3)"})
        ok = verify(s, mod=3).ok
        out.update(found_rank=s.rank, verified_mod3=ok, cyclic_invariant=is_cyclic_invariant(s),
                   orbits_used=sum(1 for o in orbits if any(e for m in o for r_ in m for e in r_)),
                   fixed_used=sum(1 for x in fixed if any(e for r_ in x for e in r_)))
        FOUND.mkdir(parents=True, exist_ok=True)
        path = FOUND / f"005_{n}x{n}x{n}_rank{s.rank}_{name}.json"
        s.save(path)
        out["scheme_file"] = str(path.relative_to(ROOT))
    if status == "UNSAT" and check_cap:
        proof = WORK / f"{name}.drat"
        t1 = time.time()
        r2 = subprocess.run(solver_cmd(solver, seed, cnf, proof), capture_output=True, text=True,
                            timeout=cap * 2)
        out["proof_solve_status"] = {10: "SAT", 20: "UNSAT"}.get(r2.returncode, f"exit{r2.returncode}")
        out["proof_bytes"] = proof.stat().st_size if proof.exists() else 0
        try:
            r3 = subprocess.run([str(DRAT_TRIM), str(cnf), str(proof), "-t", str(int(check_cap))],
                                capture_output=True, text=True, timeout=check_cap + 30)
            out["drat_trim"] = "VERIFIED" if "s VERIFIED" in r3.stdout else (
                [l for l in r3.stdout.splitlines() if l.startswith("s ")] or ["no verdict"])[0]
        except subprocess.TimeoutExpired:
            out["drat_trim"] = "TIMEOUT"
        out["proof_seconds"] = round(time.time() - t1, 2)
        if out["drat_trim"] == "VERIFIED":
            proof.unlink()  # large; reproducible from the CNF and seed
    return out


def record(res, stage, cap):
    fieldname = "GF(3)"
    best = res.get("found_rank")
    record_run(ROOT / "results" / "runs.jsonl", attempt=ATTEMPT, stage=stage,
               format=f"{res['n']}x{res['n']}x{res['n']}", field=fieldname,
               method=f"SAT {res['solver']} on cyclic-invariant Brent equations mod 3, split a={res['a']} b={res['b']}",
               best_rank_found=best if best is not None else None,
               verified=bool(res.get("verified_mod3")) if res["status"] == "SAT" else (res.get("drat_trim") == "VERIFIED"),
               budget={"max_seconds": cap, "processes": 2, "threads_each": 1}, **{"result": res})


def run_many(jobs, stage, cap, check_cap):
    results = []
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(run_instance, *j, cap, check_cap) for j in jobs]
        for f in futs:
            res = f.result()
            record(res, stage, cap)
            print(json.dumps(res), flush=True)
            results.append(res)
    return results


if __name__ == "__main__":
    stage = sys.argv[1]
    if stage == "0":
        jobs = [(2, 2, 1, "cadical", 50001), (2, 2, 1, "kissat", 50001),
                (2, 1, 3, "cadical", 50001), (2, 0, 6, "cadical", 50001),
                (2, 1, 3, "kissat", 50001), (2, 0, 6, "kissat", 50001)]
        run_many(jobs, "0", 60, 120)
    elif stage == "A":
        run_many([(3, 7, 2, "cadical", 50002), (3, 7, 2, "kissat", 50002)], "A", 1200, None)
    elif stage == "B":
        solver = sys.argv[2] if len(sys.argv) > 2 else "kissat"
        splits = [(1, 19), (2, 16), (3, 13), (4, 10), (5, 7), (6, 4), (7, 1)]
        run_many([(3, a, b, solver, 50003) for a, b in splits], "B", 600, 1200)
    else:
        raise SystemExit(__doc__)
