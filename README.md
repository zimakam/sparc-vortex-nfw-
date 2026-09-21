# SPARC VORTEX vs NFW

Reproducible analysis comparing a vortex-based dark matter halo
profile against the standard NFW profile on 171 SPARC galaxies.

## Key results

| Metric | VORTEX | NFW |
|---|---|---|
| Success rate (chi2/dof < 3) | **83.0%** | 57.9% |
| Median chi2/dof | **1.11** | 2.40 |
| Mean chi2/dof | **1.93** | 4.42 |
| Sum BIC | **8170** | 14074 |

- **Delta BIC = +5904** in favour of VORTEX
- 95% bootstrap CI: **[+3135, +9397]**
- P(Delta BIC > 0) = **1.0000**
- Head-to-head BIC wins: **141 / 171** (82.5%)

## Models

VORTEX (curl-free profile):
    v_vortex(r) = A * sqrt((x - arctan x) / x),  x = r / r_c

NFW (standard):
    v_NFW(r) = A_h * sqrt((ln(1+x) - x/(1+x)) / x),  x = r / r_s

Both use identical baryonic modelling with Upsilon_b = 1.4 Upsilon_d,
and identical parameter count (3 free parameters).

## Data

SPARC catalogue (Lelli, McGaugh & Schombert 2016, AJ, 152, 157).
Download: http://astroweb.cwru.edu/SPARC/

Files needed in data/:
- MassModels_Lelli2016c.mrt  (rotation curves)
- SPARC_Lelli2016c.mrt       (galaxy properties, Q flags)

## Install

    pip install -r requirements.txt

## Run

    python run_sparc.py        # fit both models (~3 seconds)
    python make_figures.py     # generate 4 publication figures
    python bootstrap.py        # 95% confidence intervals

## Author

Magomed K. Ziyavutdinov (Independent researcher)
zimakam@gmail.com

## License

MIT - see LICENSE.
