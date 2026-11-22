"""Attempt 006: cyclic-invariant F3 search split into cubes by matrix class (see attempts/006-cyclic-f3-cubes.md).

A cube (n, a, b, c) is CyclicF3CNF(n, a, b, pin_u1=class_reps()[c], min_class=c): orbit slot 0's u is the
representative of class c, and (minimal class form) no orbit matrix lies in class 0 (c >= 1) or is of
rank one (c >= 2). c = None is the unpinned 005 instance. Solving, decoding, verifying, saving and the
drat-trim check work as in scripts/run_005.py. Every instance gets one line in results/runs.jsonl.

usage: python3 -I scripts/run_006.py STAGE        (STAGE = A, B, C, B2 or C2)
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
from matmul.cyclic import CyclicF3CNF, class_reps, expand, is_cyclic_invariant  # noqa: E402
from matmul.record import record_run  # noqa: E402
from matmul.verify import verify  # noqa: E402

WORK = ROOT / "scratch" / "006"
FOUND = ROOT / "results" / "f3"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
ATTEMPT = "006-cyclic-f3-cubes"


def solver_cmd(solver, seed, cnf, proof=None):
    cmd = [solver, "-q", f"--seed={seed}", str(cnf)]
    if proof:
        cmd.append(str(proof))
    return cmd


def run_instance(n, a, b, c, solver, seed, cap, check_cap=None, excl=False):
    name = f"n{n}_a{a}_b{b}_" + ("free" if c is None else f"c{c}") + ("x" if excl else "") + f"_{solver}_s{seed}"
    WORK.mkdir(parents=True, exist_ok=True)
    enc = CyclicF3CNF(n, a, b) if c is None else CyclicF3CNF(n, a, b, pin_u1=class_reps(n)[c], min_class=c,
                                                             exclude_charpoly=excl)
    cnf = WORK / f"{name}.cnf"
    cnf.write_text(enc.dimacs())
    out = dict(n=n, a=a, b=b, rank_max=3 * a + b, cube=c, exclude_charpoly=excl, solver=solver, seed=seed, cap_seconds=cap,
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
        s = expand(n, orbits, fixed, {"attempt": ATTEMPT, "instance": name, "field": "GF(3)"})
        out.update(found_rank=s.rank, verified_mod3=verify(s, mod=3).ok, cyclic_invariant=is_cyclic_invariant(s),
                   orbits_used=sum(1 for o in orbits if any(e for m in o for r_ in m for e in r_)),
                   fixed_used=sum(1 for x in fixed if any(e for r_ in x for e in r_)))
        FOUND.mkdir(parents=True, exist_ok=True)
        path = FOUND / f"006_{n}x{n}x{n}_rank{s.rank}_{name}.json"
        s.save(path)
        out["scheme_file"] = str(path.relative_to(ROOT))
    if status == "UNSAT" and check_cap:
        proof = WORK / f"{name}.drat"
        t1 = time.time()
        try:
            r2 = subprocess.run(solver_cmd(solver, seed, cnf, proof), capture_output=True, text=True,
                                timeout=cap * 2)
            out["proof_solve_status"] = {10: "SAT", 20: "UNSAT"}.get(r2.returncode, f"exit{r2.returncode}")
        except subprocess.TimeoutExpired:
            out["proof_solve_status"] = "TIMEOUT"
        out["proof_bytes"] = proof.stat().st_size if proof.exists() else 0
        if out["proof_solve_status"] == "UNSAT":
            try:
                r3 = subprocess.run([str(DRAT_TRIM), str(cnf), str(proof), "-t", str(int(check_cap))],
                                    capture_output=True, text=True, timeout=check_cap + 30)
                out["drat_trim"] = "VERIFIED" if "s VERIFIED" in r3.stdout else (
                    [l for l in r3.stdout.splitlines() if l.startswith("s ")] or ["no verdict"])[0]
            except subprocess.TimeoutExpired:
                out["drat_trim"] = "TIMEOUT"
        out["proof_seconds"] = round(time.time() - t1, 2)
        if out.get("drat_trim") == "VERIFIED":
            proof.unlink()  # large; reproducible from the CNF and seed
    return out


def record(res, stage, cap):
    best = res.get("found_rank")
    cube = "unpinned" if res["cube"] is None else f"class cube {res['cube']}"
    if res.get("exclude_charpoly"):
        cube += " with char-poly exclusion"
    record_run(ROOT / "results" / "runs.jsonl", attempt=ATTEMPT, stage=stage,
               format=f"{res['n']}x{res['n']}x{res['n']}", field="GF(3)",
               method=f"SAT {res['solver']} on cyclic-invariant Brent equations mod 3, split a={res['a']} b={res['b']}, {cube}",
               best_rank_found=best,
               verified=bool(res.get("verified_mod3")) if res["status"] == "SAT" else (res.get("drat_trim") == "VERIFIED"),
               budget={"max_seconds": cap, "processes": 2, "threads_each": 1}, result=res)


def run_many(jobs, stage, cap, check_cap, excl=False):
    results = []
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        futs = [ex.submit(run_instance, *j, cap, check_cap, excl) for j in jobs]
        for f in futs:
            res = f.result()
            record(res, stage, cap)
            print(json.dumps(res), flush=True)
            results.append(res)
    return results


if __name__ == "__main__":
    stage = sys.argv[1]
    if stage == "A":
        run_many([(3, 7, 2, 0, "cadical", 60001), (3, 7, 2, 0, "cadical", 60002),
                  (3, 7, 2, None, "cadical", 60003), (3, 7, 2, None, "cadical", 60004)], "A", 300, None)
    elif stage == "B":
        run_many([(3, 2, 16, c, "cadical", 60011) for c in range(21)], "B", 120, 600)
    elif stage == "C":
        run_many([(3, 7, 1, c, "cadical", 60012) for c in range(21)], "C", 120, 600)
    elif stage == "B2":
        run_many([(3, 2, 16, c, "cadical", 60013) for c in (10, 19)], "B2", 300, 900, excl=True)
    elif stage == "C2":
        run_many([(3, 7, 1, c, "cadical", 60014) for c in range(10, 21)], "C2", 120, 600, excl=True)
    else:
        raise SystemExit(__doc__)
