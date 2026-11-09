# HPM_Duffing_Suite

This MATLAB code reproduces every number, table and figure in the manuscript
*"Accuracy, convergence and failure of the homotopy perturbation method for the
nonlinear charged particle oscillation equation: a benchmark against the exact
elliptic solution"*.

Problem: `x'' + C x + D x^3 = 0`, `x(0) = A`, `x'(0) = 0`.
Scaled form: `X'' + X + eps X^3 = 0`, with `tau = sqrt(C) t` and `eps = D A^2 / C`.

## How to run

```matlab
>> main_run_all
```

The code runs in MATLAB R2019b or later and in GNU Octave 7 or later. The published
results were produced with GNU Octave 8.4.0. A full run takes about 4 minutes, most of
it spent on the parameter map.

Outputs:
- `results/Table1..Table6_*.csv`: the data behind manuscript Tables II–VII
- `results/summary.txt`: verification checks, scaling check, map boundaries, frequency errors and timing
- `results/results.mat`: all metrics, without the time series
- `figures/Fig1..Fig5_*.png`: the manuscript figures at 300 dpi

## Files

| File | Purpose |
|---|---|
| `main_run_all.m` | Driver script: runs every step below |
| `hpm_coefficients.m` | Exact rational coefficients of the plain HPM terms u0..u4 and the HPM-FE terms X0..X4 and w1..w4 |
| `hpm_plain.m` | Plain HPM partial sums x^(n), n <= 4, with analytic derivative |
| `hpm_fe.m` | HPM with frequency expansion (He / Lindstedt–Poincaré type), order <= 3 |
| `exact_cn.m` | Exact solution `cn(sqrt(1+eps) tau | eps/(2(1+eps)))`, plus exact frequency and period |
| `energy_nd.m` | Hamiltonian `H = X'^2/2 + X^2/2 + eps X^4/4` |
| `hpm_hierarchy_ode.m` | Numerical integration of the HPM recurrences, used for verification |
| `verify_terms.m` | Independent checks of the closed-form terms (agreement < 1e-12) |
| `case_definitions.m` | The seven (C, D) parameter sets |
| `compute_case_metrics.m` | MAE, RMSE, error ratios, energy drift and frequency error per case |
| `parameter_map.m` | 81 x 60 map over (eps, number of periods) |
| `timing_study.m` | Wall-clock cost of each method |
| `write_tables.m` | CSV output |
| `make_figures.m` | Figures 1–5 |

## Notes

- The exact solution is the reference throughout. The `ode45` solution (RelTol 1e-12,
  AbsTol 1e-14) agrees with it to within 4e-12.
- `hpm_fe` returns NaN when the truncated frequency series gives Omega^2 <= 0.
  At order 3 this happens for eps >= 4.48.
- Absolute timings depend on the software. Octave's `ode45` is much slower than MATLAB's.
