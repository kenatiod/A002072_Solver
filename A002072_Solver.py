#!/usr/bin/env python3
"""
A002072_Solver.py — version 10
=================================
Compute OEIS A002072: largest m such that m and m+1 are both p_r-smooth
(i.e., all prime factors of m*(m+1) lie in {p_1, p_2, ..., p_r}).

Method (Stormer 1897 / Lehmer 1964) and all mathematical identities are
unchanged from version 9 — see the identity block below.  Version 10 is
a systems rework that removes the v9 feeder bottleneck.

Version 10 changes over v9
---------------------------
1. IN-WORKER DISCRIMINANT GENERATION.  v9 generated every squarefree
   subset product D in a feeder process and shipped them to workers in
   batches; at r = 32+ the feeder DFS and queue traffic became the
   dominant cost (v9 wall time doubled per r — pure 2^r task overhead).
   v10 partitions the subset space by the membership pattern of the TOP
   H primes: a task is a single integer (the high mask), and each
   worker runs the pruned DFS over the LOW (smaller) primes itself,
   with the high-product as the DFS root.  IPC per task: one small int
   out, one aggregate tuple back.  Blocks whose high-product already
   exceeds the prefilter threshold are discarded in O(1) with their
   2^L subsets counted in closed form.

2. POOL ARCHITECTURE.  The hand-rolled task/result queues, feeder
   process, sentinels, and shared progress counters are replaced by a
   multiprocessing.Pool with a per-process initializer (which performs
   the cypari2/gp setup once per worker) and imap_unordered over the
   block indices.  Fewer moving parts, identical semantics.

3. ACCOUNTING preserved exactly:
       sum(block n_d) + sum(block prefiltered) == 2^r - 1
   is checked and reported as "subset_accounting_ok".  Prefilter
   threshold remains the corrected (2*max_m + 1)^2 of v9.

4. CLI: --batch_size is replaced by --block_log2 (default 20): each
   block is ~2^block_log2 subsets.  JSON: "n_batches"/"batch_size" are
   replaced by "n_blocks"/"block_log2"/"split_H".

Key identities (unchanged from v9; proofs in the v9 docstring)
---------------------------------------------------------------
With x^2 - D*y^2 = 1 and D a product of primes of P_r (smooth by
construction):
  (I1) m+1 = x^2 smooth  <=>  x smooth;  omega(m+1) = omega(x).
  (I2) m = D*y^2 smooth  <=>  y smooth;  omega(m) = |primes(D) U primes(y)|.
  (I3) x odd: m = (x-1)/2, m+1 = (x+1)/2 are BOTH smooth <=> y smooth.
  (I4) y1 | yn (Lucas divisibility) => if y1 is not smooth, the
       discriminant yields no candidate at any index: skip the loop.
Iteration depth L2 = 2*max(3, 2r+4, (p_r+1)/2) subsumes the removed
4*D equations (even-y solutions form an index <= 2 subgroup).

Semantics preserved from v8/v9
-------------------------------
* Candidates are recorded even when m > max_m (the cap bounds the
  SEARCH, not the recording).
* "n_smooth_pairs_found" counts discriminants that produced at least
  one smooth pair.
* Prime-complete detection: omega(m) + omega(m+1) == r.
* Adaptive max_m resolution, A002072_KNOWN override, JSON audit trail,
  stderr capture, cypari2/subprocess-gp dual path: unchanged.

Usage examples
--------------
# Unlimited run for r=1..15 (verify against OEIS):
python A002072_Solver_10.py --start_r 1 --end_r 15

# Adaptive run for r=16..35 (expo_margin=6 recommended):
python A002072_Solver_10.py --start_r 16 --end_r 35 --expo_margin 6

# Fixed cap for a specific r:
python A002072_Solver_10.py --start_r 28 --end_r 28 --max_m_expo 29

By Ken Clements, June 2026.
Based on the Lehmer/Stormer approach coded by Lucas A. Brown (stormer.py).
Version 9/10 changes co-developed with Claude (Anthropic).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from multiprocessing import get_context, cpu_count
from typing import Dict, List, Optional, Tuple

try:
    sys.set_int_max_str_digits(0)
except AttributeError:
    pass

PROGRAM_NAME    = "A002072_Solver"
PROGRAM_VERSION = 10

# ---------------------------------------------------------------------------
# OEIS A002072 ground truth (r=1..27, unconditionally verified)
# ---------------------------------------------------------------------------
A002072_KNOWN: Dict[int, int] = {
    1:  1,
    2:  8,
    3:  80,
    4:  4374,
    5:  9800,
    6:  123200,
    7:  336140,
    8:  11859210,
    9:  11859210,
    10: 177182720,
    11: 1611308699,
    12: 3463199999,
    13: 63927525375,
    14: 421138799639,
    15: 1109496723125,
    16: 1453579866024,
    17: 20628591204480,
    18: 31887350832896,
    19: 31887350832896,
    20: 119089041053696,
    21: 2286831727304144,
    22: 9591468737351909375,
    23: 9591468737351909375,
    24: 9591468737351909375,
    25: 9591468737351909375,
    26: 9591468737351909375,
    27: 19316158377073923834000,
}

# ---------------------------------------------------------------------------
# Optional fast-path: gmpy2 and cypari2
# ---------------------------------------------------------------------------
try:
    import gmpy2
    from gmpy2 import mpz as _mpz
    def mpz(x): return _mpz(x)
    HAS_GMPY2 = True
except ImportError:
    def mpz(x): return int(x)  # type: ignore[misc]
    HAS_GMPY2 = False

CYPARI2_SESSION = None  # set per-worker in the pool initializer

try:
    import cypari2 as _cypari2
    HAS_CYPARI2 = True
except ImportError:
    HAS_CYPARI2 = False

try:
    _popcount = int.bit_count          # unbound method: _popcount(mask)
    _popcount(0)
except (AttributeError, TypeError):
    def _popcount(x: int) -> int:      # type: ignore[misc]
        return bin(x).count("1")

# ---------------------------------------------------------------------------
# Small utilities
# ---------------------------------------------------------------------------

def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()

def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()

def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

def primes_up_to_r(r: int) -> List[int]:
    primes: List[int] = []
    c = 2
    while len(primes) < r:
        if all(c % p for p in primes):
            primes.append(c)
        c += 1 if c == 2 else 2
    return primes

def env_block(script_path: str, argv: List[str], gp_path: str) -> Dict:
    src = open(script_path, "rb").read()
    block: Dict = {
        "script_path":   os.path.abspath(script_path),
        "script_sha256": sha256_bytes(src),
        "program":       f"{PROGRAM_NAME} v{PROGRAM_VERSION}",
        "command_line":  " ".join(argv),
        "python_version": sys.version.replace("\n", " "),
        "platform":      platform.platform(),
        "gmpy2":         str(HAS_GMPY2),
        "cypari2":       str(HAS_CYPARI2),
        "gp_path":       gp_path,
    }
    try:
        r = subprocess.run([gp_path, "-q"], input="print(version());quit\n",
                           capture_output=True, text=True, timeout=10)
        block["gp_version"] = (r.stdout or "").strip().splitlines()[-1]
    except Exception as e:
        block["gp_version"] = f"unavailable ({e})"
    try:
        r = subprocess.run(["git", "rev-parse", "HEAD"],
                           capture_output=True, text=True)
        block["git_commit"] = r.stdout.strip()
    except Exception:
        block["git_commit"] = None
    return block

# ---------------------------------------------------------------------------
# Adaptive max_m resolution (unchanged from v9)
# ---------------------------------------------------------------------------

def resolve_max_m(r: int, json_dir: str, expo_margin: int,
                  max_m_expo: int) -> Tuple[int, int, str]:
    """
    Returns (max_m, effective_expo, source_description).
    max_m == 0 means unlimited.
    """
    if max_m_expo > 0:
        return (10 ** max_m_expo, max_m_expo, "fixed")

    if expo_margin == 0:
        if r in A002072_KNOWN:
            known = A002072_KNOWN[r]
            expo  = int(math.log10(known)) + 2
            print(f"  [unlimited->capped] r={r}: A002072_KNOWN[{r}]={known:,} "
                  f"-> max_m = 10^{expo}")
            return (10 ** expo, expo, f"A002072_KNOWN[{r}] unlimited-override")
        return (0, 0, "unlimited")

    prev = r - 1
    if prev in A002072_KNOWN:
        anchor = A002072_KNOWN[prev]
        source = f"A002072_KNOWN[{prev}]"
    else:
        jf = os.path.join(json_dir, f"smooth_r{prev:02d}.json")
        if not os.path.isfile(jf):
            print(f"  [adaptive] No prior JSON for r={prev}; running unlimited.")
            return (0, 0, f"unlimited (no prior JSON for r={prev})")
        with open(jf) as f:
            d = json.load(f)
        if not d.get("search_complete", False):
            print(f"  [adaptive] Prior JSON r={prev} is NOT complete; "
                  f"running unlimited for safety.")
            return (0, 0, f"unlimited (r={prev} search was incomplete)")
        anchor = d["last_smooth_m"]
        source = f"smooth_r{prev:02d}.json"

    if anchor < 1:
        return (0, 0, f"unlimited (anchor={anchor} is invalid)")

    expo  = int(math.log10(anchor)) + expo_margin + 1
    max_m = 10 ** expo
    print(f"  [adaptive] r={r}: anchor={anchor:,} [{source}] "
          f"-> max_m = 10^{expo}  (floor(log10)+{expo_margin}+1)")
    return (max_m, expo, source)

# ---------------------------------------------------------------------------
# PARI / Pell interface (cypari2 preferred, subprocess gp fallback)
# ---------------------------------------------------------------------------

_PELL_GP_SRC = r"""
pellxy(D, max_x=0)={
  if(D<=0, error("D<=0"));
  if(issquare(D), error("D is square"));
  my(a0=sqrtint(D), m=0, d=1, a=a0, p0=1, p1=a0, q0=0, q1=1);
  while(p1^2 - D*q1^2 != 1,
    m = d*a - m;
    d = (D - m^2)/d;
    a = (a0 + m)\d;
    my(p2=a*p1+p0, q2=a*q1+q0);
    p0=p1; p1=p2; q0=q1; q1=q2;
    if(max_x>0 && p1>max_x, return([0,0]));
  );
  if(max_x>0 && p1>max_x, return([0,0]));
  [p1, q1];
};
"""

_GP_PROC:         Optional[subprocess.Popen] = None
_GP_PATH_GLOBAL:  str  = "gp"
_VEC2_RE = re.compile(r"^\[\s*(\d+)\s*,\s*(\d+)\s*\]\s*$")
_SENTINEL = "__A002072_SOLVER_SEP__"

def _gp_kill():
    global _GP_PROC
    if _GP_PROC is None:
        return
    try:
        _GP_PROC.kill()
    except Exception:
        pass
    _GP_PROC = None

def _gp_start() -> subprocess.Popen:
    p = subprocess.Popen(
        [_GP_PATH_GLOBAL, "-q"],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, bufsize=1,
    )
    assert p.stdin and p.stdout
    p.stdin.write(_PELL_GP_SRC + "\n")
    p.stdin.write(f'print("{_SENTINEL}_start");\n')
    p.stdin.write("print(pellxy(46));\n")
    p.stdin.write("v=pellxy(46); print(v[1]^2 - 46*v[2]^2);\n")
    p.stdin.write("print(pellxy(46, 100));\n")
    p.stdin.write(f'print("{_SENTINEL}_end");\n')
    p.stdin.flush()
    lines, in_block = [], False
    while True:
        line = p.stdout.readline()
        if not line:
            raise RuntimeError("gp handshake: EOF")
        s = line.strip()
        if s == f"{_SENTINEL}_start":
            in_block = True
            continue
        if s == f"{_SENTINEL}_end":
            break
        if in_block:
            lines.append(s)
    if len(lines) < 3 or lines[1] != "1" or lines[2] not in ("[0, 0]", "[0,0]"):
        raise RuntimeError(f"gp handshake failed: {lines}")
    return p

def _gp_eval(expr: str) -> str:
    global _GP_PROC
    if _GP_PROC is None:
        _GP_PROC = _gp_start()
    assert _GP_PROC.stdin and _GP_PROC.stdout
    _GP_PROC.stdin.write(f'print("{_SENTINEL}_start");\n')
    _GP_PROC.stdin.write(f"print({expr});\n")
    _GP_PROC.stdin.write(f'print("{_SENTINEL}_end");\n')
    _GP_PROC.stdin.flush()
    lines, in_block = [], False
    while True:
        line = _GP_PROC.stdout.readline()
        if not line:
            raise RuntimeError("gp: EOF (process died)")
        s = line.strip()
        if s == f"{_SENTINEL}_start":
            in_block = True
            continue
        if s == f"{_SENTINEL}_end":
            break
        if in_block:
            lines.append(s)
    return "\n".join(lines).strip()

def _pell_fundamental(D: int, max_x: int = 0,
                      retries: int = 3) -> Tuple[int, int]:
    """
    Return (x1, y1), the fundamental solution to x^2-Dy^2=1, or (0, 0)
    if x1 > max_x (when max_x > 0).  cypari2 fast path with subprocess
    gp fallback.  The Pell identity is verified by the caller in mpz.
    """
    if HAS_CYPARI2 and CYPARI2_SESSION is not None:
        pari = CYPARI2_SESSION
        try:
            if max_x > 0:
                v = pari(f"pellxy({D}, {max_x})")
            else:
                v = pari(f"pellxy({D})")
            return int(v[0]), int(v[1])
        except Exception:
            pass  # fall through to subprocess

    last_exc: Optional[Exception] = None
    for _ in range(retries):
        try:
            expr = f"pellxy({D}, {max_x})" if max_x > 0 else f"pellxy({D})"
            out  = _gp_eval(expr)
            if "***" in out or "error" in out.lower():
                raise RuntimeError(f"gp error: {out[:200]}")
            m = _VEC2_RE.match(out)
            if not m:
                raise RuntimeError(f"Unexpected gp output: {out!r}")
            return int(m.group(1)), int(m.group(2))
        except Exception as e:
            last_exc = e
            _gp_kill()
    raise RuntimeError(f"Pell solver failed for D={D}: {last_exc}")

# ---------------------------------------------------------------------------
# Worker: per-process state, initializer, and the block task
# ---------------------------------------------------------------------------

_S: Dict = {}   # per-process solver state, installed by _solver_init

def _solver_init(primes: List[int], gp_path: str, max_m: int,
                 L2: int, H: int, verbose: bool) -> None:
    """
    Pool initializer: runs once per worker process.  Sets up the
    cypari2 session (silently) or the subprocess-gp fallback, and
    installs all per-r constants in the module-level dict _S.
    """
    global _GP_PATH_GLOBAL, CYPARI2_SESSION
    try:
        sys.set_int_max_str_digits(0)
    except AttributeError:
        pass
    _GP_PATH_GLOBAL = gp_path
    if HAS_CYPARI2:
        try:
            pari = _cypari2.Pari()
            # cypari2's allocatemem() does a Python-level print() of
            # "PARI stack size set to ..." on sys.stdout (NOT the C
            # stderr fd).  silent=True suppresses it; fall back to
            # capturing stdout for older cypari2 versions.
            try:
                pari.allocatemem(256 * 1024 * 1024, silent=True)
            except TypeError:
                import io, contextlib
                with contextlib.redirect_stdout(io.StringIO()):
                    pari.allocatemem(256 * 1024 * 1024)
            pari(_PELL_GP_SRC)
            v = pari("pellxy(46)")
            assert int(v[0]) ** 2 - 46 * int(v[1]) ** 2 == 1
            CYPARI2_SESSION = pari
        except Exception as e:
            CYPARI2_SESSION = None
            print(f"  [worker] cypari2 init failed: {e}; using subprocess gp",
                  flush=True)

    r = len(primes)
    L = r - H
    primes_mpz = [mpz(p) for p in primes]

    if HAS_GMPY2:
        _remove = gmpy2.remove
    else:
        def _remove(n, p):  # type: ignore[misc]
            k = 0
            while n % p == 0:
                n //= p
                k += 1
            return n, k

    ONE = mpz(1)

    def factor_smooth(n) -> Tuple[bool, int]:
        """Trial-divide n (mpz >= 1) by P_r -> (is_smooth, dividing-mask)."""
        mask = 0
        for i in range(r):
            p = primes_mpz[i]
            if n % p == 0:
                mask |= 1 << i
                n, _ = _remove(n, p)
                if n <= ONE:
                    return True, mask
        return n == ONE, mask

    _S.update({
        "primes":        primes,
        "r":             r,
        "H":             H,
        "L":             L,
        "primes_low":    primes[:L],
        "primes_high":   primes[L:],
        "max_m":         max_m,
        "max_x":         (2 * max_m + 1) if max_m > 0 else 0,
        "T":             (2 * max_m + 1) ** 2 if max_m > 0 else 0,
        "L2":            L2,
        "ONE":           ONE,
        "factor_smooth": factor_smooth,
        "verbose":       verbose,
    })

def _process_D(D: int, agg: Dict) -> None:
    """
    Run the full v9 per-discriminant algorithm on one D, mutating the
    block aggregate `agg`.  (Pell solve -> identity check -> y1 gate ->
    L2-deep recurrence with y/x-based candidate extraction.)
    """
    r             = _S["r"]
    max_x         = _S["max_x"]
    L2            = _S["L2"]
    ONE           = _S["ONE"]
    factor_smooth = _S["factor_smooth"]

    x1, y1 = _pell_fundamental(D, max_x=max_x)
    if x1 == 0 and y1 == 0:
        agg["n_cut"] += 1
        return

    Dm  = mpz(D)
    x1m = mpz(x1)
    y1m = mpz(y1)
    if x1m * x1m - Dm * y1m * y1m != ONE:
        raise RuntimeError(f"Pell identity check failed for D={D}")

    # y1 gate (I4)
    sm_y1, _unused = factor_smooth(y1m)
    if not sm_y1:
        agg["n_gate"] += 1
        return

    sm_D, mask_D = factor_smooth(Dm)
    if not sm_D:
        raise RuntimeError(f"D={D} not smooth (generator bug)")

    Dy1m = Dm * y1m
    x, y = x1m, y1m
    found_pair = False

    for _ in range(L2):
        if max_x and x > max_x:
            break
        agg["n_iter"] += 1

        sm_y, mask_y = factor_smooth(y)
        if sm_y:
            # ---- half candidate: m = (x-1)/2 (x odd) -- both smooth by I3
            if x & 1:
                mm = (x - 1) >> 1
                sm_a, mask_a = factor_smooth(mm)
                sm_b, mask_b = factor_smooth(mm + 1)
                if not (sm_a and sm_b):
                    raise RuntimeError(
                        f"identity I3 violated at D={D} (internal error)")
                found_pair = True
                if mm > agg["best_m"]:
                    agg["best_m"] = mm
                if _popcount(mask_a) + _popcount(mask_b) == r:
                    agg["pc_hits"].append(int(mm))
            # ---- primary candidate: m = D*y^2, m+1 = x^2 (I1, I2)
            sm_x, mask_x = factor_smooth(x)
            if sm_x:
                mm = Dm * y * y          # formed only on a hit
                found_pair = True
                if mm > agg["best_m"]:
                    agg["best_m"] = mm
                if _popcount(mask_D | mask_y) + _popcount(mask_x) == r:
                    agg["pc_hits"].append(int(mm))

        x, y = x1m * x + Dy1m * y, x1m * y + y1m * x

    if found_pair:
        agg["n_pair_d"] += 1

def _solver_block(high_mask: int):
    """
    One task: the fixed high-prime membership pattern `high_mask` over
    the top H primes.  The worker computes the base product, then runs
    the pruned DFS over the LOW primes (ascending, incremental
    products, prune when the partial product exceeds T = (2*max_m+1)^2,
    counting the pruned subtree as 2^(L-j) - 1 in closed form).

    Block accounting (checked globally by the collector):
      n_d + n_prefiltered == 2^L      (high_mask > 0)
      n_d + n_prefiltered == 2^L - 1  (high_mask == 0; the global empty
                                       subset D=1 is excluded)

    Returns:
      (high_mask, n_d, n_prefiltered, n_cut, n_gate, n_iter,
       best_m_int, n_pair_d, pc_hits, n_err, err_samples)
    """
    try:
        primes_low  = _S["primes_low"]
        primes_high = _S["primes_high"]
        L           = _S["L"]
        T           = _S["T"]

        base = 1
        for i in range(len(primes_high)):
            if (high_mask >> i) & 1:
                base *= primes_high[i]

        agg: Dict = {
            "n_d": 0, "pref": 0, "n_cut": 0, "n_gate": 0, "n_iter": 0,
            "best_m": mpz(0), "n_pair_d": 0, "pc_hits": [],
            "n_err": 0, "err_samples": [],
        }

        def run_D(D: int) -> None:
            agg["n_d"] += 1
            try:
                _process_D(D, agg)
            except Exception as e:
                agg["n_err"] += 1
                if len(agg["err_samples"]) < 5:
                    agg["err_samples"].append(f"D={D}: {e!r}")

        if high_mask > 0 and T and base > T:
            # Entire block exceeds the prefilter threshold: every low
            # subset (including the empty one, D = base) is pruned.
            agg["pref"] = 1 << L
        else:
            if high_mask > 0:
                run_D(base)            # the empty-low subset, D = base
            def rec(i: int, prod: int) -> None:
                for j in range(i, L):
                    np_ = prod * primes_low[j]
                    if T and np_ > T:
                        agg["pref"] += (1 << (L - j)) - 1
                        break
                    run_D(np_)
                    rec(j + 1, np_)
            rec(0, base)

        return (high_mask, agg["n_d"], agg["pref"], agg["n_cut"],
                agg["n_gate"], agg["n_iter"], int(agg["best_m"]),
                agg["n_pair_d"], agg["pc_hits"], agg["n_err"],
                agg["err_samples"])

    except Exception as e:
        return (high_mask, 0, 0, 0, 0, 0, 0, 0, [], 1,
                [f"block {high_mask}: {e!r}"])

# ---------------------------------------------------------------------------
# run_one_r  (v10: Pool over high-mask blocks)
# ---------------------------------------------------------------------------

def run_one_r(
    r:              int,
    primes:         List[int],
    max_m:          int,
    effective_expo: int,
    anchor_source:  str,
    anchor_value:   int,
    expo_margin:    int,
    n_workers:      int,
    block_log2:     int,
    json_dir:       str,
    gp_path:        str,
    verbose:        bool,
) -> Dict:
    p_r       = primes[-1]
    start_utc = utc_now()
    t0        = time.time()

    L_base = max(3, 2 * r + 4, (p_r + 1) // 2)
    L2     = 2 * L_base

    # Block split: H high bits per task, ~2^block_log2 subsets per task,
    # with a floor guaranteeing at least ~8 tasks per worker.
    if r <= 1:
        H = 0
    else:
        H = max(r - block_log2,
                math.ceil(math.log2(max(2, n_workers * 8))))
        H = max(0, min(H, r - 1))
    num_blocks = 1 << H

    ctx = get_context("fork") if sys.platform == "darwin" else get_context()

    last_smooth_m  = 0
    n_smooth_pairs = 0
    n_errors       = 0
    n_d_run        = 0
    n_prefiltered  = 0
    n_cutoffs      = 0
    n_gate_skips   = 0
    n_solutions    = 0
    n_collected    = 0
    error_samples: List[str] = []
    all_pc_hits:   set       = set()

    last_progress_time = t0

    with ctx.Pool(processes=n_workers,
                  initializer=_solver_init,
                  initargs=(primes, gp_path, max_m, L2, H, verbose)) as pool:
        for res in pool.imap_unordered(_solver_block, range(num_blocks),
                                       chunksize=1):
            (hm, n_d, pref, n_cut, n_gate, n_it,
             b_m, n_pair, pchits, n_e, errsamp) = res
            n_collected   += 1
            n_d_run       += n_d
            n_prefiltered += pref
            n_cutoffs     += n_cut
            n_gate_skips  += n_gate
            n_solutions   += n_it
            n_smooth_pairs += n_pair
            n_errors      += n_e
            if b_m > last_smooth_m:
                last_smooth_m = b_m
            if pchits:
                all_pc_hits.update(pchits)
            for s in errsamp:
                if len(error_samples) < 20:
                    error_samples.append(s)

            now = time.time()
            if verbose and (now - last_progress_time) >= 30.0:
                print(f"  blocks={n_collected:,}/{num_blocks:,}"
                      f"  D_run={n_d_run:,}"
                      f"  prefiltered={n_prefiltered:,}"
                      f"  gate_skips={n_gate_skips:,}"
                      f"  best_m={last_smooth_m:,}"
                      f"  pc_hits={len(all_pc_hits)}"
                      f"  elapsed={(now - t0)/60:.1f}min",
                      flush=True)
                last_progress_time = now

    n_subsets     = (1 << r) - 1
    accounting_ok = (n_d_run + n_prefiltered == n_subsets) \
                    and (n_collected == num_blocks)
    if not accounting_ok:
        print(f"  [!] r={r}: subset accounting mismatch: "
              f"run={n_d_run:,} + prefiltered={n_prefiltered:,} "
              f"!= {n_subsets:,} (blocks {n_collected}/{num_blocks})")

    if verbose or n_prefiltered > 0:
        print(f"  r={r} p_r={p_r} discriminants_run={n_d_run:,} "
              f"(+{n_prefiltered:,} pre-filtered by D>(2*max_m+1)^2) "
              f"cutoffs={n_cutoffs:,} gate_skips={n_gate_skips:,} "
              f"max_m={'unlimited' if max_m == 0 else f'10^{effective_expo}'} "
              f"workers={n_workers} blocks={num_blocks} L2={L2}")

    runtime         = time.time() - t0
    end_utc         = utc_now()
    search_complete = (n_errors == 0) and accounting_ok
    pc_hits_sorted  = sorted(all_pc_hits)
    n_pc_hits       = len(pc_hits_sorted)

    known_value  = A002072_KNOWN.get(r)
    xcheck_known = None
    if known_value is not None:
        xcheck_known = (last_smooth_m == known_value)
        if not xcheck_known:
            print(f"  [MISMATCH] r={r}: computed {last_smooth_m:,} "
                  f"but A002072_KNOWN[{r}]={known_value:,}")
        elif verbose:
            print(f"  [OK] r={r}: matches A002072_KNOWN[{r}]={known_value:,}")

    if verbose and n_pc_hits > 0:
        print(f"  [PC] r={r}: {n_pc_hits} prime-complete hit(s): "
              + ", ".join(f"{m:,}" for m in pc_hits_sorted))
    elif verbose:
        print(f"  [PC] r={r}: no prime-complete pairs found")

    summary = {
        "program":             f"{PROGRAM_NAME} v{PROGRAM_VERSION}",
        "r":                   r,
        "p_r":                 p_r,
        "primes":              primes,
        "last_smooth_m":       last_smooth_m,
        "last_smooth_m_plus1": last_smooth_m + 1 if last_smooth_m > 0 else 0,
        "effective_max_m_expo": effective_expo,
        "expo_margin":         expo_margin,
        "anchor_source":       anchor_source,
        "anchor_value":        anchor_value,
        "search_complete":     search_complete,
        "n_discriminants":     n_d_run,
        "n_prefiltered_subsets": n_prefiltered,
        "n_subsets_total":     n_subsets,
        "subset_accounting_ok": accounting_ok,
        "n_pell_cutoffs":      n_cutoffs,
        "n_y1_gate_skips":     n_gate_skips,
        "n_solutions_tried":   n_solutions,
        "n_smooth_pairs_found": n_smooth_pairs,
        "n_discriminant_errors": n_errors,
        "error_samples":       error_samples,
        "prime_complete_hits": pc_hits_sorted,
        "n_prime_complete_hits": n_pc_hits,
        "lehmer_L_base":       L_base,
        "l2_iteration_depth":  L2,
        "block_log2":          block_log2,
        "split_H":             H,
        "n_blocks":            num_blocks,
        "runtime_seconds":     round(runtime, 3),
        "start_utc":           start_utc,
        "end_utc":             end_utc,
        "a002072_known_value": known_value,
        "xcheck_known_matches": xcheck_known,
        "workers":             n_workers,
    }

    out_path = os.path.join(json_dir, f"smooth_r{r:02d}.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    if verbose:
        print(f"  -> {out_path}  last_smooth_m={last_smooth_m:,}"
              f"  runtime={runtime/60:.2f} min"
              f"  complete={search_complete}"
              f"  pc_hits={n_pc_hits}")

    return summary

# ---------------------------------------------------------------------------
# Stderr capture (unchanged from v9)
# ---------------------------------------------------------------------------

def _redirect_stderr_to_tmpfile() -> Tuple[int, str]:
    import tempfile
    stderr_fd = sys.stderr.fileno()
    saved_fd  = os.dup(stderr_fd)
    tmp_fd, tmp_path = tempfile.mkstemp(prefix="a002072_stderr_", suffix=".txt")
    os.dup2(tmp_fd, stderr_fd)
    os.close(tmp_fd)
    return saved_fd, tmp_path

def _restore_stderr(saved_fd: int, tmp_path: str,
                    verbose: bool = False) -> None:
    stderr_fd = sys.stderr.fileno()
    os.dup2(saved_fd, stderr_fd)
    os.close(saved_fd)
    if not os.path.isfile(tmp_path):
        return
    size = os.path.getsize(tmp_path)
    if size == 0:
        os.remove(tmp_path)
        return
    with open(tmp_path, "r", errors="replace") as f:
        content = f.read().strip()
    os.remove(tmp_path)
    if content:
        print(f"\n[stderr captured during run]\n{content}",
              file=sys.stderr, flush=True)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    ap = argparse.ArgumentParser(
        description=f"{PROGRAM_NAME} v{PROGRAM_VERSION}: extend OEIS A002072"
    )
    ap.add_argument("--version", action="store_true",
                    help="Print version info as JSON and exit.")
    ap.add_argument("--start_r", type=int, default=1)
    ap.add_argument("--end_r",   type=int, default=27)
    ap.add_argument("--max_m_expo", type=int, default=0,
                    help="Fixed cap: max_m = 10^max_m_expo.  0 = unlimited.")
    ap.add_argument("--expo_margin", type=int, default=0,
                    help="Adaptive cap exponent margin (recommended: 6).")
    ap.add_argument("--workers",  type=int, default=0,
                    help="Worker processes (default: all CPUs).")
    ap.add_argument("--block_log2", type=int, default=20,
                    help="Target log2(subsets) per worker block (default 20).")
    ap.add_argument("--outdir",   default="A002072_Solver_runs")
    ap.add_argument("--gp_path",  default="gp")
    ap.add_argument("--skip_done", action="store_true")
    ap.add_argument("--verbose",  action="store_true")
    args = ap.parse_args()

    if args.version:
        info = env_block(__file__, sys.argv, args.gp_path)
        info["a002072_known_r_max"] = max(A002072_KNOWN)
        print(json.dumps(info, indent=2, sort_keys=True))
        return

    if args.max_m_expo > 0 and args.expo_margin > 0:
        print(f"[warning] Both --max_m_expo={args.max_m_expo} and "
              f"--expo_margin={args.expo_margin} set.  "
              f"Using fixed --max_m_expo; ignoring --expo_margin.")

    n_workers = args.workers if args.workers > 0 else (cpu_count() or 1)
    ensure_dir(args.outdir)

    try:
        subprocess.run([args.gp_path, "-q"], input="print(1);quit\n",
                       text=True, capture_output=True, check=True, timeout=15)
    except Exception as e:
        print(f"[!] Cannot run gp at '{args.gp_path}': {e}")
        print("    Install PARI/GP or pass --gp_path /path/to/gp")
        sys.exit(1)

    _stderr_saved_fd: int = 0
    _stderr_tmp_path: str = ""
    if not args.verbose:
        _stderr_saved_fd, _stderr_tmp_path = _redirect_stderr_to_tmpfile()

    import atexit
    def _atexit_restore():
        if _stderr_saved_fd:
            _restore_stderr(_stderr_saved_fd, _stderr_tmp_path,
                            verbose=args.verbose)
    atexit.register(_atexit_restore)

    print(f"[+] {PROGRAM_NAME} v{PROGRAM_VERSION}")
    print(f"[+] gmpy2: {HAS_GMPY2}   cypari2: {HAS_CYPARI2}")
    print(f"[+] workers={n_workers}  block_log2={args.block_log2}  "
          f"outdir={args.outdir}")
    cap_mode = ("unlimited"
                if args.max_m_expo == 0 and args.expo_margin == 0
                else (f"fixed 10^{args.max_m_expo}" if args.max_m_expo > 0
                      else f"adaptive (expo_margin={args.expo_margin})"))
    print(f"[+] cap mode: {cap_mode}")
    print(f"[+] r range: {args.start_r}..{args.end_r}")
    print(f"[+] v10: in-worker block generation + pool architecture "
          f"(v9 math unchanged)")
    print()

    env      = env_block(__file__, sys.argv, args.gp_path)
    env_path = os.path.join(args.outdir, "environment.json")
    with open(env_path, "w") as f:
        json.dump(env, f, indent=2, sort_keys=True)

    for r in range(args.start_r, args.end_r + 1):
        out_path = os.path.join(args.outdir, f"smooth_r{r:02d}.json")
        if args.skip_done and os.path.isfile(out_path):
            print(f"[=] r={r}: already done ({out_path}), skipping.")
            continue

        primes = primes_up_to_r(r)

        max_m, effective_expo, anchor_source = resolve_max_m(
            r           = r,
            json_dir    = args.outdir,
            expo_margin = args.expo_margin,
            max_m_expo  = args.max_m_expo,
        )
        anchor_value = A002072_KNOWN.get(r - 1, 0)

        print(f"[>] r={r}  p_r={primes[-1]}  "
              f"{'unlimited' if max_m == 0 else f'max_m=10^{effective_expo}'}")

        summary = run_one_r(
            r              = r,
            primes         = primes,
            max_m          = max_m,
            effective_expo = effective_expo,
            anchor_source  = anchor_source,
            anchor_value   = anchor_value,
            expo_margin    = args.expo_margin,
            n_workers      = n_workers,
            block_log2     = args.block_log2,
            json_dir       = args.outdir,
            gp_path        = args.gp_path,
            verbose        = args.verbose,
        )

        status = ("MATCHES OEIS"  if summary.get("xcheck_known_matches")
                  else ("NEW RESULT" if summary.get("a002072_known_value") is None
                        else "!!! MISMATCH !!!"))
        print(f"    last_smooth_m = {summary['last_smooth_m']:,}")
        print(f"    runtime       = {summary['runtime_seconds']/60:.2f} min  "
              f"complete={summary['search_complete']}  {status}")
        print()

    print(f"[+] Done.  Results in {args.outdir}/")


if __name__ == "__main__":
    if sys.platform == "darwin":
        import multiprocessing as _mp
        try:
            _mp.set_start_method("fork")
        except RuntimeError:
            pass
    main()
