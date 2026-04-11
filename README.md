Contrast-Space Projection (CSP) for Network Meta-Analysis
This repository provides the official R implementation of the Contrast-Space Projection (CSP) framework, an exact and invariant study-based decomposition of direct and indirect evidence contributions in Network Meta-Analysis (NMA).

Overview
The CSP framework addresses the "reproducibility gap" in NMA. While traditional methods often struggle to reproduce NMA estimates from marginal pairwise summaries—particularly when multi-arm trials induce complex within-study correlations—CSP provides an exact, invariant, and study-based decomposition.

Key Features
Exact Decomposition: Reproduces NMA estimates precisely using a projection-based approach.

Invariance: Results are independent of the choice of baseline treatments or contrast definitions.

Multi-Arm Integration: Naturally handles within-study correlations without data loss.

Visual Diagnostics: Includes interactive tools for assessing network "tension" and evidence contribution.

Repository Structure
csp_functions.R: The core library containing the mathematical engine for the CSP framework.

csp_example.R: A walkthrough script using an illustrative dataset (Treatments A–E, Studies I–V) to demonstrate the analysis pipeline.

Getting Started
Prerequisites
Ensure you have R installed along with the following packages:

R
install.packages(c("Matrix", "MASS", "dplyr", "ggplot2", "plotly"))
Running the Example
Place csp_functions.R and csp_example.R in the same directory.

Open R or RStudio and run:

R
source("csp_example.R")
The script will output the Projection Matrix (P), the Global Q Inconsistency Test, and generate interactive visualizations for evidence decomposition.

Citation
If you use this code or the CSP framework in your research, please cite both the software repository and the original manuscript:

Article
Wang, C., Zhang, Y., Jin, Z., & O'Connor, A. (2026). Contrast-Space Projection for Network Meta-Analysis: An Exact and Invariant Study-Based Decomposition of Direct and Indirect Contributions. Research Synthesis Methods.

Software
Wang, C., Zhang, Y., Jin, Z., & O'Connor, A. (2026). CSP-NMA: R implementation of the Contrast-Space Projection framework. https://github.com/your-username/CSP-NMA.

BibTeX
Code snippet
@article{wang2026contrast,
  title={Contrast-Space Projection for Network Meta-Analysis: An Exact and Invariant Study-Based Decomposition of Direct and Indirect Contributions},
  author={Wang, Chong and Zhang, Yanqi and Jin, Zhezhen and O'Connor, Annette},
  journal={Research Synthesis Methods},
  year={2026}
}

@misc{wang2026cspcode,
  author={Wang, Chong and Zhang, Yanqi and Jin, Zhezhen and O'Connor, Annette},
  title={CSP-NMA: R implementation of the Contrast-Space Projection framework},
  year={2026},
  publisher={GitHub},
  journal={GitHub repository},
  howpublished={\url{https://https://github.com/your-username/CSP-NMA}}
}
License
This project is licensed under the MIT License.

Copyright (c) 2026 Chong Wang, Yanqi Zhang, Zhezhen Jin, and Annette O'Connor.
