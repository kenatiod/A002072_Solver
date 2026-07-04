# A002072 Solver

Certified computational extension tools for [OEIS A002072](https://oeis.org/A002072): the largest integer $m$ such that $m$ and $m+1$ are both $p_r$-smooth, where $p_r$ is the $r$-th prime.

Equivalently, for each order $r$, this program searches for the largest consecutive smooth pair

m, m+1 with $P^+(m(m+1)) <= p_r$.

The solver is part of a broader project on **prime-complete products of two consecutive integers**, where a product $m(m+1)$ is prime-complete of order $r$ when

$rad(m(m+1)) = p_r\#$.

The A002072 values provide the smooth ceiling $L_r$; if a prime-complete product of order $r$ exists, then necessarily $m <= L_r$.

---

## Current status

The repository contains the solver and run archives under:

```text
A002072_Solver_runs/
```

The latest run archive extends A002072 beyond the previously known OEIS range. In particular, the current database includes the certified value

```text
A002072(37) = 27129807647978258459761875
```

for $p_37 = 157$.

Selected recent values:

| r | $p_r$ | A002072(r) |
|---:|---:|---:|
| 27 | 103 | 19316158377073923834000 |
| 28 | 107 | 19316158377073923834000 |
| 29 | 109 | 19316158377073923834000 |
| 30 | 113 | 19316158377073923834000 |
| 31 | 127 | 19316158377073923834000 |
| 32 | 131 | 19316158377073923834000 |
| 33 | 137 | 124225935845233319439173 |
| 34 | 139 | 124225935845233319439173 |
| 35 | 149 | 124225935845233319439173 |
| 36 | 151 | 124225935845233319439173 |
| 37 | 157 | 27129807647978258459761875 |

The $r=37$ result is especially useful for the prime-complete project because

$L_37(L_37+1)$

jumps across multiple primorial intervals while still remaining far below $P_37 = p_37\#$.

---

## Mathematical method

The solver uses the classical Størmer--Lehmer reduction from consecutive smooth integers to Pell equations.

For squarefree $D$ supported on the first $r$ primes, solve

```text
$x^2 - D y^2 = 1$.
```

The relevant identities are:

```text
$m     = D y^2$,
$m + 1 = x^2$,
```

and, for odd `$x$`, the half-solution branch

```text
$m     = (x - 1)/2$,
$m + 1 = (x + 1)/2$.
```

The implementation iterates Pell solutions on each squarefree discriminant $D | p_r\#$, checks smoothness against the first $r$ primes, records the largest smooth pair found, and separately records prime-complete hits when the distinct prime support has size exactly $r$.

The code also uses the Lucas divisibility gate:

```text
$y_1 | y_n$.
```

Therefore, if the fundamental $y_1$ is not $p_r$-smooth, the entire Pell branch can be skipped.

---

## Version 10 overview

`A002072_Solver.py` version 10 keeps the mathematics of version 9 unchanged and reworks the task architecture for large $r$.

The major change is **in-worker discriminant generation**. Instead of generating every squarefree subset product $D$ in a feeder process and shipping large batches to workers, version 10 partitions the subset space by the membership pattern of the top $H$ primes. Each worker receives a high-prime mask and performs the pruned low-prime DFS locally.

This reduces interprocess communication from many large discriminant batches to a small integer task plus one aggregate result.

Version 10 also replaces the earlier hand-rolled queue architecture with `multiprocessing.Pool`, preserves exact subset accounting, and writes a JSON audit trail for every completed $r$.

---

## Requirements

Required:

- Python 3.10+ recommended
- [PARI/GP](https://pari.math.u-bordeaux.fr/) available as `gp`

Optional but recommended:

- `gmpy2` for faster big integer arithmetic
- `cypari2` for faster in-process PARI calls

The program automatically falls back to subprocess `gp` when `cypari2` is unavailable.

On macOS with Homebrew, a typical setup is:

```bash
brew install pari
python3 -m pip install gmpy2
# cypari2 is optional and may require PARI headers/libraries configured correctly
python3 -m pip install cypari2
```

Check your environment with:

```bash
python3 A002072_Solver.py --version
```

---

## Quick start

Verify small known values:

```bash
python3 A002072_Solver.py --start_r 1 --end_r 15
```

Adaptive extension run, using prior completed JSON outputs as anchors:

```bash
python3 A002072_Solver.py --start_r 16 --end_r 37 --expo_margin 6
```

Run a single fixed-cap calculation:

```bash
python3 A002072_Solver.py --start_r 37 --end_r 37 --max_m_expo 30
```

Use all CPUs by default, or specify worker count:

```bash
python3 A002072_Solver.py --start_r 37 --end_r 37 --expo_margin 6 --workers 10
```

Skip previously completed JSON outputs:

```bash
python3 A002072_Solver.py --start_r 16 --end_r 37 --expo_margin 6 --skip_done
```

Choose the output directory:

```bash
python3 A002072_Solver.py --start_r 34 --end_r 37 --expo_margin 6 --outdir A002072_Solver_runs
```

---

## Important command-line options

| Option | Meaning |
|---|---|
| `--start_r` | First order `r` to run. |
| `--end_r` | Last order `r` to run. |
| `--max_m_expo` | Fixed cap `max_m = 10^max_m_expo`; `0` means no fixed cap. |
| `--expo_margin` | Adaptive cap margin. Recommended value for extension runs: `6`. |
| `--workers` | Number of worker processes. Default: all CPUs. |
| `--block_log2` | Target log2 subset count per high-mask block. Default: `20`. |
| `--outdir` | Directory for JSON outputs and `environment.json`. |
| `--gp_path` | Path to `gp` if it is not on `PATH`. |
| `--skip_done` | Skip `smooth_rXX.json` files that already exist. |
| `--verbose` | Print progress during block collection. |
| `--version` | Print environment/version metadata as JSON and exit. |

---

## Adaptive cap mode

With `--expo_margin M`, the solver uses the previous completed value as an anchor and sets

```text
$max_m = 10^(floor(log10(anchor)) + M + 1)$.
```

For example, with `--expo_margin 6`, each new run searches several orders of magnitude beyond the previous certified A002072 value.

Important: `max_m` bounds the **search**, not the recording. If a candidate is encountered beyond `max_m` before the Pell branch is cut off, it is still recorded. The cap is used to make the Pell/discriminant search finite and auditable.

If the previous JSON is missing or incomplete, the solver falls back to an unlimited run for safety.

---

## Output files

Each completed order writes a JSON file:

```text
A002072_Solver_runs/smooth_rXX.json
```

Each JSON output contains:

- $r$, $p_r$, and the prime list
- `last_smooth_m`
- `last_smooth_m_plus1`
- adaptive cap metadata
- `search_complete`
- subset/discriminant accounting
- Pell cutoffs
- $y_1$ gate skips
- number of Pell solution iterations tried
- number of smooth-pair discriminants found
- prime-complete hits, if any
- runtime and UTC timestamps
- known-value cross-checks when applicable

The run directory also contains:

```text
environment.json
```

which records the script hash, command line, Python version, platform, optional library availability, PARI/GP version, and git commit when available.

---

## Certification checks

For every order $r$, the solver checks exact subset accounting:

```text
n_discriminants + n_prefiltered_subsets == $2^r - 1$.
```

The JSON field

```text
subset_accounting_ok
```

must be `true` for a run to be considered complete.

A run is considered complete only when:

```text
search_complete == true
```

which requires zero discriminant errors and valid subset accounting.

For orders already known to OEIS, the solver cross-checks the computed value against the internal `A002072_KNOWN` table and reports whether it matches.

---

## Prime-complete hit detection

In addition to computing A002072 values, the solver checks for prime-complete products.

A smooth pair is recorded as prime-complete when the number of distinct prime factors across $m$ and $m+1$ is exactly $r$:

```text
omega(m) + omega(m+1) == r.
```

Prime-complete hits are written to the JSON field:

```text
prime_complete_hits
```

with count:

```text
n_prime_complete_hits
```

For the broader prime-complete project, these fields provide an independent audit of whether any smooth-pair search also found a prime-complete product in the certified range.

---

## Repository layout

Suggested repository layout:

```text
A002072_Solver/
  A002072_Solver.py
  README.md
  A002072_Solver_runs/
    environment.json
    smooth_r16.json
    smooth_r17.json
    ...
    smooth_r37.json
```

If publishing or citing a run archive, prefer a tagged release or commit hash rather than the moving `main` branch.

Recommended citation form:

```text
Ken Clements, A002072_Solver, release <tag>, commit <hash>,
https://github.com/kenatiod/A002072_Solver
```

---

## Reproducibility notes

For a referee-facing or archival release:

1. Tag the repository.
2. Record the commit hash.
3. Preserve all `smooth_rXX.json` files.
4. Preserve `environment.json`.
5. Generate SHA256 hashes for all run outputs.
6. Record the exact command lines used.
7. Archive the release externally if possible, for example with Zenodo.

A simple hash manifest can be generated with:

```bash
find A002072_Solver_runs -type f -print0 | sort -z | xargs -0 shasum -a 256 > SHA256SUMS.txt
```

---

## Relationship to prime-complete consecutive products

Let

```text
$L_r$ = A002072(r).
```

If $m(m+1)$ is prime-complete of order $r$, then both $m$ and $m+1$ are $p_r$-smooth, so

```text
$m <= L_r$.
```

Thus A002072 provides a smooth ceiling for possible prime-complete products.

The larger prime-complete project combines this ceiling with other finite-certification and tail-elimination methods, including:

- direct value-axis searches,
- Lehmer--Clements/Pell branch elimination,
- CRT pruning,
- RAD-CRT floor comparisons,
- and structural reduction theorems for the remaining tail.

This repository supplies the A002072 smooth-ceiling layer and records any prime-complete hits encountered during that search.

---

## Acknowledgments

The mathematical method follows the classical work of Størmer and Lehmer on consecutive smooth integers and Pell equations.

This implementation is by Ken Clements, June 2026.

It is based in part on the Lehmer/Størmer approach coded by Lucas A. Brown (`stormer.py`). Version 9/10 systems changes were co-developed with Claude (Anthropic), while the present repository is maintained by Ken Clements.

