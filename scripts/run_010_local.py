"""Attempt 010, stream "local": S3 rank-23 pool by local search, and rank-22 local non-existence theorems.

usage:
  python3 -I scripts/run_010_local.py lns ITER CAP CPU_LIMIT      local search on the pool (seeds 130001 + i)
  python3 -I scripts/run_010_local.py thm KEY K CAP CPU_LIMIT     rank-22 local theorem: pool scheme KEY (class key
                                                                   prefix), every K-subset of its S3 slots frozen
  python3 -I scripts/run_010_local.py one KEY '[frozen slots]' '[a,b,c,d]' CAP SEED   re-run one thm instance
Every solver instance gets one line in results/runs.jsonl (stage "L"); stdout is a JSON-lines log.

lns: pick a pool scheme and 3 or 4 of its S3 slots at random, free them (same slot types), forbid the freed
orbits from all coming back, solve at the same rank (cap CAP). SAT answers are verified mod 3 and checked for
S3 invariance; a new class key is added to the pool and saved as results/f3/010_local_*.json. A SAT answer of
rank < 23 stops everything.

thm: for every K-subset F of the scheme's S3 slots with rank(F) <= 22, and every S3 split (a, b, c, d) of
22 - rank(F) (every slot may be zero, so the exact-rank splits cover all smaller ones), solve
S3LocalF3CNF(frozen = F, free = (a, b, c, d)). UNSAT is re-solved with a DRAT proof and checked by drat-trim.
All UNSAT and VERIFIED for every K-subset means: no S3-invariant scheme of rank <= 22 over F3 contains K of
the scheme's S3 orbits.
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import random
import shutil
import sys
import time
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from matmul.invariants import class_key  # noqa: E402
from matmul.record import record_run  # noqa: E402
from matmul.s3 import is_s3_invariant, scheme_s3_split  # noqa: E402
from matmul.s3_local import S3LocalF3CNF, group_slots, s3_slots, slot_list, slots_rank  # noqa: E402
from matmul.satrun import solve_instance  # noqa: E402
from matmul.scheme import Scheme  # noqa: E402

WORK = ROOT / "scratch" / "010" / "local"
FOUND = WORK / "found"
POOL = WORK / "pool.jsonl"
DRAT_TRIM = ROOT / "vendor" / "drat-trim" / "drat-trim"
RUNS = ROOT / "results" / "runs.jsonl"
ATTEMPT = "010-s3-f3"
S3_006 = ROOT / "results" / "f3" / "010_3x3x3_rank23_s3_from006.json"


def load_pool():
    pool = [dict(class_key="eebc48753638175a", file=str(S3_006.relative_to(ROOT)), how="006 scheme in S3 form")]
    if POOL.exists():
        pool += [json.loads(l) for l in POOL.read_text().splitlines() if l.strip()]
    return pool


def split_of(items):
    return tuple(len(x) for x in group_slots(items))


def splits_of_rank(r):
    if r == 0:
        return [(0, 0, 0, 0)]
    out = []
    for a in range(r // 6 + 1):
        for b in range((r - 6 * a) // 3 + 1):
            for c in range((r - 6 * a - 3 * b) // 2 + 1):
                out.append((a, b, c, r - 6 * a - 3 * b - 2 * c))
    return out


def run_one(enc, name, seed, cap, check_cap, meta):
    t = time.time()
    res = solve_instance(enc, name, "cadical", seed, cap, work=WORK, found=FOUND, root=ROOT, drat_trim=DRAT_TRIM,
                         attempt=ATTEMPT, prefix="010_local", check_cap=check_cap, meta=meta)
    res["wall"] = round(time.time() - t, 2)
    cnf = WORK / f"{name}.cnf"
    if cnf.exists() and (res["status"] != "UNSAT" or res.get("drat_trim") == "VERIFIED"):
        cnf.unlink()  # reproducible from the script, the frozen slots and the seed
    if res["status"] == "SAT":
        s = Scheme.load(ROOT / res["scheme_file"])
        res.update(s3_invariant=is_s3_invariant(s), class_key=class_key(s, 3))
        if res["s3_invariant"]:
            res["s3_split_found"] = list(scheme_s3_split(s))
    return res


def record(res, cap, method):
    sat = res["status"] == "SAT"
    record_run(RUNS, attempt=ATTEMPT, stage="L", format="3x3x3", field="GF(3)", method=method,
               best_rank_found=res.get("found_rank"),
               verified=bool(res.get("verified_mod3") and res.get("s3_invariant")) if sat else res.get("drat_trim") == "VERIFIED",
               budget={"max_seconds": cap, "processes": 2, "threads_each": 1}, result=res)


def cpu(res):
    return res.get("seconds", 0) + res.get("proof_seconds", 0)


def lns(iters, cap, cpu_limit):
    rng = random.Random(130001)
    pool = load_pool()
    keys = {p["class_key"] for p in pool}
    used = 0.0

    def job(i):
        p = rng.choice(pool)
        items = slot_list(s3_slots(Scheme.load(ROOT / p["file"])))
        m = rng.choice((3, 3, 4))
        idx = sorted(rng.sample(range(len(items)), m))
        frozen = [items[j] for j in range(len(items)) if j not in idx]
        free = [items[j] for j in idx]
        enc = S3LocalF3CNF(3, group_slots(frozen), split_of(free), forbid=group_slots(free))
        seed = 130001 + i
        name = f"lns_{p['class_key'][:8]}_free{''.join(map(str, idx))}_s{seed}"
        meta = dict(mode="lns", source=p["class_key"], free_slots=idx, free_split=list(split_of(free)),
                    frozen_split=list(split_of(frozen)), s3_rank_max=23)
        return run_one(enc, name, seed, cap, 300, meta)

    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        i = 0
        futs = set()
        while i < iters or futs:
            while len(futs) < 2 and i < iters and used < cpu_limit:
                futs.add(ex.submit(job, i))
                i += 1
            if not futs:
                break
            done, futs = cf.wait(futs, return_when=cf.FIRST_COMPLETED)
            for f in done:
                res = f.result()
                used += cpu(res)
                new = False
                if res["status"] == "SAT" and res.get("verified_mod3") and res.get("s3_invariant"):
                    if res["found_rank"] < 23:
                        res["STOP"] = "rank < 23 found"
                    if res["class_key"] not in keys:
                        new = True
                        keys.add(res["class_key"])
                        dst = ROOT / "results" / "f3" / f"010_local_3x3x3_rank{res['found_rank']}_{res['class_key']}.json"
                        shutil.copy(ROOT / res["scheme_file"], dst)
                        entry = dict(class_key=res["class_key"], file=str(dst.relative_to(ROOT)),
                                     how=f"lns from {res['source']} freeing slots {res['free_slots']} seed {res['seed']}",
                                     s3_split=res["s3_split_found"], rank=res["found_rank"])
                        pool.append(entry)
                        with POOL.open("a") as fh:
                            fh.write(json.dumps(entry) + "\n")
                res["new_class"] = new
                res["cpu_used"] = round(used, 1)
                record(res, cap, f"SAT cadical, S3 local search: S3 rank-23 scheme {res['source']} with slots "
                                 f"{res['free_slots']} free (same types), freed orbits forbidden from all returning")
                print(json.dumps(res), flush=True)
                if res.get("STOP"):
                    raise SystemExit("rank < 23 found")


def thm_one(key, frozen_idx, sp, cap, seed, check_cap=1200):
    """Re-run one local-theorem instance (for example a timeout) with another cap and seed."""
    p = next(q for q in load_pool() if q["class_key"].startswith(key))
    items = slot_list(s3_slots(Scheme.load(ROOT / p["file"])))
    F = tuple(frozen_idx)
    frozen = [items[j] for j in F]
    assert slots_rank(frozen) + 6 * sp[0] + 3 * sp[1] + 2 * sp[2] + sp[3] == 22
    enc = S3LocalF3CNF(3, group_slots(frozen), tuple(sp))
    name = f"thm_{p['class_key'][:8]}_k{len(F)}_keep{''.join(map(str, F))}_free{''.join(map(str, sp))}_s{seed}"
    meta = dict(mode="thm", retry=True, source=p["class_key"], k=len(F), frozen_slots=list(F),
                frozen_rank=slots_rank(frozen), free_split=list(sp), s3_rank_max=22)
    res = run_one(enc, name, seed, cap, check_cap, meta)
    record(res, cap, f"SAT cadical, S3 rank-22 local theorem (retry): S3 slots {list(F)} of {p['class_key']} "
                     f"frozen, free S3 split {tuple(sp)}")
    print(json.dumps(res), flush=True)


def thm(key, k, cap, cpu_limit, check_cap=600):
    pool = load_pool()
    p = next(q for q in pool if q["class_key"].startswith(key))
    s = Scheme.load(ROOT / p["file"])
    items = slot_list(s3_slots(s))
    jobs = []
    for F in combinations(range(len(items)), k):
        frozen = [items[j] for j in F]
        r = 22 - slots_rank(frozen)
        if r < 0:
            continue
        for sp in splits_of_rank(r):
            jobs.append((F, sp))
    print(json.dumps(dict(mode="thm-plan", source=p["class_key"], k=k, slots=[ty for ty, _ in items],
                          subsets=len({F for F, _ in jobs}), instances=len(jobs))), flush=True)
    used = 0.0

    def job(F, sp):
        frozen = [items[j] for j in F]
        enc = S3LocalF3CNF(3, group_slots(frozen), sp)
        seed = 131001
        name = f"thm_{p['class_key'][:8]}_k{k}_keep{''.join(map(str, F))}_free{''.join(map(str, sp))}_s{seed}"
        meta = dict(mode="thm", source=p["class_key"], k=k, frozen_slots=list(F), frozen_rank=slots_rank(frozen),
                    free_split=list(sp), s3_rank_max=22)
        return run_one(enc, name, seed, cap, check_cap, meta)

    tally = {}
    with cf.ThreadPoolExecutor(max_workers=2) as ex:
        it = iter(jobs)
        futs = set()
        while True:
            while len(futs) < 2 and used < cpu_limit:
                nxt = next(it, None)
                if nxt is None:
                    break
                futs.add(ex.submit(job, *nxt))
            if not futs:
                break
            done, futs = cf.wait(futs, return_when=cf.FIRST_COMPLETED)
            for f in done:
                res = f.result()
                used += cpu(res)
                res["cpu_used"] = round(used, 1)
                tag = res["status"] + ("/" + str(res.get("drat_trim")) if res["status"] == "UNSAT" else "")
                tally[tag] = tally.get(tag, 0) + 1
                record(res, cap, f"SAT cadical, S3 rank-22 local theorem: S3 slots {res['frozen_slots']} of {res['source']} "
                                 f"frozen, free S3 split {tuple(res['free_split'])}")
                print(json.dumps(res), flush=True)
                if res["status"] == "SAT" and res.get("verified_mod3") and res.get("s3_invariant") and res["found_rank"] <= 22:
                    dst = ROOT / "results" / "f3" / f"010_local_3x3x3_rank{res['found_rank']}_{res['class_key']}.json"
                    shutil.copy(ROOT / res["scheme_file"], dst)
                    print(json.dumps(dict(STOP="rank-22 S3 scheme found", file=str(dst))), flush=True)
                    raise SystemExit("rank 22 found")
    print(json.dumps(dict(mode="thm-summary", source=p["class_key"], k=k, instances=len(jobs), tally=tally,
                          cpu=round(used, 1), complete=sum(tally.values()) == len(jobs))), flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "lns":
        lns(int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]))
    elif mode == "thm":
        thm(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
    elif mode == "one":
        thm_one(sys.argv[2], json.loads(sys.argv[3]), json.loads(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]))
    else:
        raise SystemExit(__doc__)
