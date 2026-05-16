# CLP Output Comparison: Old vs New (pilot87)

## Old CLP (pre-facelift)

```
Coin LP version devel, build May 15 2026
command line - clp-old /home/hgambinisantos/inst/netlib/pilot87.mps.gz
Model was imported from /home/hgambinisantos/inst/netlib/pilot87.mps.gz in 0.007024 seconds
Perturbing problem by 0.001% of 0.23783515 - largest nonzero change 0.00010007932 ( 3955.6806%) - largest zero change 0.00010006998
3315  Obj 160.12866 Primal inf 16925074 (881)
5089  Obj 281.52188 Primal inf 11501476 (816)
6588  Obj 303.39119 Primal inf 206985.51 (624)
Optimal - objective value 301.71065
Optimal objective 301.7106548 - 7324 iterations time 2.502, Presolve 0.01
```

## New CLP (deterministic output: `-progress 0 -progressIter 1000`)

```
CLP devel (git:ae7cbb15) — COIN-OR LP Solver
  args: /home/hgambinisantos/inst/netlib/pilot87.mps.gz -progress 0 -progressIter 1000 -solve

Problem: PILOT87
  |Rows| =  2030   |Cols| =  4883   |NZ| = 73152

  Coefficient ranges:
    Matrix    ∈ [1e-06,   1000]   κ ≈ 1e+09
    Cost      ∈ [1e-06, 3.2963]   κ ≈ 3.296e+06
    Bounds    ∈ [ 0.02,  16475]   κ ≈ 8.237e+05
    RHS       ∈ [ 0.01,   6500]   κ ≈ 650000

  Perturbation: 0.001% of 0.237835
  CLP presolve: 2030 rows, 4883 cols → 1784 rows, 4406 cols (0.01s)

▶ LP solve

  ──────── ──────── ─────────────── ──────────── ──────────── ────────
  Phase        Iter       Objective   Primal inf     Dual inf     Time
  ──────── ──────── ─────────────── ──────────── ──────────── ────────
  Dual            1   -2.700448e+03   4.3965e+08   3.3760e+00    0.019
  Dual         1001    4.081575e+01   8.5848e+07   0.0000e+00    0.115
  Dual         2001    1.233961e+02   7.7717e+06   0.0000e+00    0.324
  Dual         3001    1.508299e+02   8.1721e+07   0.0000e+00    0.615
  Dual         4001    2.190579e+02   1.6739e+07   0.0000e+00    0.977
  Dual         5001    2.769474e+02   1.5379e+08   0.0000e+00     1.43
  Dual         6001    2.993924e+02   1.5982e+05   0.0000e+00     1.90
  Dual         7001    3.044514e+02   8.5984e+04   0.0000e+00     2.35
  Primal       7369    3.017108e+02   0.0000e+00   0.0000e+00     2.55
  ──────── ──────── ─────────────── ──────────── ──────────── ────────

✔ Optimal — Obj: 301.711   Iters: 7369   Time: 2.55s
```

## New CLP (default, time-based progress every ~0.7s)

```
CLP devel (git:ae7cbb15) — COIN-OR LP Solver
  args: /home/hgambinisantos/inst/netlib/pilot87.mps.gz

Problem: PILOT87
  |Rows| =  2030   |Cols| =  4883   |NZ| = 73152

  Coefficient ranges:
    Matrix    ∈ [1e-06,   1000]   κ ≈ 1e+09
    Cost      ∈ [1e-06, 3.2963]   κ ≈ 3.296e+06
    Bounds    ∈ [ 0.02,  16475]   κ ≈ 8.237e+05
    RHS       ∈ [ 0.01,   6500]   κ ≈ 650000

  Perturbation: 0.001% of 0.237835
  CLP presolve: 2030 rows, 4883 cols → 1784 rows, 4406 cols (0.01s)

▶ LP solve

  ──────── ──────── ─────────────── ──────────── ──────────── ────────
  Phase        Iter       Objective   Primal inf     Dual inf     Time
  ──────── ──────── ─────────────── ──────────── ──────────── ────────
  Dual            1   -2.700448e+03   4.3965e+08   3.3760e+00    0.019
  Dual         3238    1.579738e+02   1.8666e+07   0.0000e+00    0.719
  Dual         4980    2.742172e+02   1.5379e+08   0.0000e+00     1.43
  Dual         6523    3.030264e+02   2.2898e+05   0.0000e+00     2.13
  Primal       7369    3.017108e+02   0.0000e+00   0.0000e+00     2.54
  ──────── ──────── ─────────────── ──────────── ──────────── ────────

✔ Optimal — Obj: 301.711   Iters: 7369   Time: 2.54s
```

## Key Differences

| Feature | Old | New |
|---------|-----|-----|
| Problem summary | Import time only | Dimensions, coefficient ranges, condition numbers |
| Presolve info | Only time in final line | Full: original → presolved dimensions + time |
| Perturbation | Full verbose message | Concise one-liner |
| Progress format | Unstructured lines | Aligned table with Phase/Iter/Obj/Inf/Time |
| Progress frequency | Time-based (~0.7s) | Time-based (default) or iteration-based (`-progressIter N`) |
| Iteration 0 | Not shown | Shows first actual iteration |
| Final status | Two lines (status + timing) | Single summary line with ✔ |
| Deterministic mode | Not available | `-progress 0 -progressIter N` |
