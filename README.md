# Contrast-Space Projection (CSP) for Network Meta-Analysis

GitHub repository: https://github.com/ChongWangStat/CSP-NMA

This repository provides the official R implementation of the **Contrast-Space Projection (CSP)** framework, an exact and invariant study-based decomposition of direct and indirect evidence contributions in Network Meta-Analysis (NMA).

---

## Overview

The CSP framework addresses the "reproducibility gap" in NMA. While traditional methods often struggle to reproduce NMA estimates from marginal pairwise summaries—particularly when multi-arm trials induce complex within-study correlations—CSP provides an exact, invariant, and study-based decomposition.

---

## Key Features

* **Exact Decomposition:** Reproduces NMA estimates precisely using a projection-based approach.  
* **Invariance:** Results are independent of the choice of baseline treatments or contrast definitions.  
* **Multi-Arm Integration:** Naturally handles within-study correlations without data loss.  
* **Visual Diagnostics:** Includes tools for assessing network "tension" and evidence contribution.  
* **Flexible Model Specification:** Supports both fixed-effects and random-effects NMA within the same projection framework.

---

## Model Specification

The implementation supports both:

- **Fixed-effects NMA** (default): uses within-study covariance only  
- **Random-effects NMA**: incorporates a heterogeneity parameter \( \tau^2 \) into the covariance structure  

Under random-effects models:
- \( \tau^2 \) is estimated (e.g., via REML)
- The covariance matrix \( \mathbf V \) is constructed using \( \hat{\tau}^2 \)
- All subsequent projection and CSP decomposition procedures remain unchanged conditional on \( \mathbf V \)

---

## Repository Structure

* **`csp_functions.R`**: Core library implementing CSP estimation and decomposition  
* **`csp_example.R`**: Walkthrough script using an illustrative dataset (Treatments A–E, Studies S1–S5)

---

## Getting Started

### Prerequisites

Install required R packages:

```r
install.packages(c("Matrix", "MASS", "dplyr", "ggplot2", "plotly"))