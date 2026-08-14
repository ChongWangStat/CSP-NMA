# Contrast-Space Projection (CSP) for Network Meta-Analysis

If you use this code, please cite the associated manuscript.

GitHub repository: https://github.com/ChongWangStat/CSP-NMA

This repository provides the official R implementation of the Contrast-Space Projection (CSP) framework, an exact and invariant study-based decomposition of direct and indirect evidence contributions in Network Meta-Analysis (NMA).

## Overview

The CSP framework addresses the reproducibility gap in NMA. The fitted NMA treatment effects are obtained from the same generalized least-squares model used for the analysis; CSP does not create a second set of network estimates. Instead, for every target comparison, the canonical direct, indirect, study-level, and path-level contributions reconstruct the fitted NMA estimate exactly, up to numerical precision.

For a fixed-effects NMA, the implementation also provides the standard generalized Cochran Q decomposition:

`Q_net = Q_het + Q_inc`,

where `Q_net` is the global network Q statistic, `Q_het` measures within-design heterogeneity, and `Q_inc` measures between-design inconsistency. These are complementary diagnostics from the same NMA analysis rather than competing sets of results.

For a random-effects NMA, heterogeneity is modeled directly through `tau^2`. Conditional on the fitted covariance matrix `V(tau^2)`, the corresponding algebraic diagnostic decomposition is reported as

`Q_net_RE = Q_error + Q_inc_RE`,

where `Q_error` is residual within-design error after accounting for the modeled heterogeneity and `Q_inc_RE` is the remaining between-design inconsistency component conditional on the fitted `tau^2`. These random-effects quantities are descriptive conditional diagnostics; the fixed-effects chi-square p-values are not reported for them.

## Key Features

- **Exact decomposition:** Reproduces each fitted NMA estimate from its canonical direct and indirect study/path contributions.
- **Invariance:** Results are independent of equivalent within-study contrast parameterizations.
- **Multi-arm integration:** Retains the covariance structure induced by multi-arm studies.
- **Canonical study/path representation:** Removes within-study algebraic redundancy before assembling indirect network paths.
- **Generalized Cochran Q diagnostics:** For fixed effects, reports `Q_net`, `Q_het`, and `Q_inc` and checks `Q_net = Q_het + Q_inc`. For random effects, reports `tau2`, `Q_net_RE`, `Q_error`, and `Q_inc_RE` and checks `Q_net_RE = Q_error + Q_inc_RE`.
- **Reproducibility checks:** Verifies that the projection, aggregated direct/indirect contributions, and full path decomposition reproduce the same fitted NMA estimate.
- **Visual diagnostics:** Includes forest, tension, path-based, and three-dimensional contribution displays.
- **Flexible modeling:** Supports fixed-effects and random-effects NMA within the same projection framework.

## Model Specification

The implementation supports both fixed-effects and random-effects NMA. Under the common-heterogeneity random-effects model, `tau^2` is estimated by REML unless supplied by the user, and each multi-arm heterogeneity block has contrast variance `tau^2` with the corresponding shared-arm covariance structure. The REML search interval is expanded automatically if the optimum lies near its initial upper boundary. Subsequent NMA estimation and CSP contribution decompositions remain algebraically exact conditional on the fitted covariance matrix `V(tau^2)`.

The interpretation of the Q decomposition differs by model. Under fixed effects, `Q_het` is the standard within-design heterogeneity statistic and the chi-square tests for `Q_net`, `Q_het`, and `Q_inc` are reported. Under random effects, heterogeneity has already been modeled through `tau^2`; therefore the same within-design quadratic is renamed `Q_error` and interpreted as residual error after accounting for modeled heterogeneity. The random-effects output is `tau2`, `Q_net_RE`, `Q_error`, and `Q_inc_RE`, satisfying the exact conditional identity `Q_net_RE = Q_error + Q_inc_RE`. Ordinary fixed-effects chi-square p-values are intentionally not reported for these random-effects components.

For formal inconsistency inference under random effects, a separate random-effects design-by-treatment interaction model should be fitted and the added inconsistency parameters tested jointly (for example, with a global Wald test). That separate inferential model is not currently implemented in this repository.

## Repository Structure

- `csp_functions.R`: Core mathematical and graphical implementation of CSP, including the generalized Cochran Q decomposition and exact reproducibility checks.
- `csp_example.R`: Walkthrough using an illustrative dataset (Treatments A-E, Studies S1-S5).

## Getting Started

### Prerequisites

Install R and the required packages:

```r
install.packages(c("Matrix", "MASS", "dplyr", "ggplot2", "plotly"))
```

### Running the Example

Place `csp_functions.R` and `csp_example.R` in the same directory, open R or RStudio, and run:

```r
source("csp_example.R")
```

The script runs the fixed-effects illustrative analysis by default and outputs the projection matrix, the generalized Cochran Q decomposition (`Q_net`, `Q_het`, and `Q_inc`), exact reproducibility checks, and the illustrative CSP decomposition and visualizations. A commented random-effects block is included as an optional example; users may uncomment it to obtain `tau2`, `Q_net_RE`, `Q_error`, and `Q_inc_RE`.

For the included fixed-effects example, the Q diagnostics are approximately:

- `Q_net = 5.6868`, `df = 5`, `p = 0.3379`
- `Q_het = 0.3570`, `df = 1`, `p = 0.5502`
- `Q_inc = 5.3298`, `df = 4`, `p = 0.2551`

and `Q_net - Q_het - Q_inc` is zero up to numerical precision.

If the optional random-effects block in `csp_example.R` is uncommented, the results are approximately:

- `tau2 = 0.0034683`
- `Q_net_RE = 4.5932`
- `Q_error = 0.3416`
- `Q_inc_RE = 4.2516`

and `Q_net_RE - Q_error - Q_inc_RE` is zero up to numerical precision. Here `tau2` is the modeled heterogeneity variance; `Q_error` is not a second heterogeneity test.

## Citation

If you use this code or the CSP framework in your research, please cite both the software repository and the original manuscript:

### Article

Wang, C., Zhang, Y., Jin, Z., & O'Connor, A. (2026). *Contrast-Space Projection for Network Meta-Analysis: An Exact and Invariant Study-Based Decomposition of Direct and Indirect Contributions.*

### Software

Wang, C., Zhang, Y., Jin, Z., & O'Connor, A. (2026). *CSP-NMA: R implementation of the Contrast-Space Projection framework.* GitHub repository.

## BibTeX

```bibtex
@article{wang2026contrast,
  title={Contrast-Space Projection for Network Meta-Analysis: An Exact and Invariant Study-Based Decomposition of Direct and Indirect Contributions},
  author={Wang, Chong and Zhang, Yanqi and Jin, Zhezhen and O'Connor, Annette},
  journal={},
  year={2026}
}

@misc{wang2026cspcode,
  author={Wang, Chong and Zhang, Yanqi and Jin, Zhezhen and O'Connor, Annette},
  title={CSP-NMA: R implementation of the Contrast-Space Projection framework},
  year={2026},
  publisher={GitHub},
  journal={GitHub repository},
  howpublished={\url{https://github.com/ChongWangStat/CSP-NMA}}
}
```

## License

This project is licensed under the MIT License.

Copyright (c) 2026 Chong Wang, Yanqi Zhang, Zhezhen Jin, and Annette O'Connor.
