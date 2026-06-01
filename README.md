# S-CARD-CMSA for the IEEE CEC 2026 Multimodal Optimization Competition

This repository contains the MATLAB source code, method description, and supporting materials for **S-CARD-CMSA**, a score-aware candidate-archive and density-filtered reporting framework developed for the **IEEE CEC 2026 Competition on Benchmarking Niching Methods for Multimodal Optimization**.

S-CARD-CMSA is built on **RS-CMSA-ESII** and introduces two conservative extensions:

1. **Passive secondary candidate archive**  
   Stores restart-level best candidates without changing the original RS-CMSA-ESII search process.

2. **Score-aware density-filtered final reporting**  
   Constructs the final reported solution set by balancing robust peak ratio (RPR) and precision-driven F1-score.

The method preserves the original RS-CMSA-ESII sampling, covariance adaptation, taboo-region update, restart logic, hill-valley checking, and termination mechanisms. The proposed changes are applied only to candidate retention and final reporting.

---

## Repository Contents

A typical repository structure is:

```text
.
├── README.md
├── data/
|   ├── fstar-data.csv
|   ├── num-normal.csv
|   ├── num-uniform.csv
|   ├── pid-data.xlsx
|   ├── sequences.csv
├── src/
│   ├── S_CARD_CSMA.m
│   ├── main.m
│   └── utility functions
├── TR/
│   └── CEC_2026_Technical_Report.pdf
├── results/
│   ├── submission_log.csv
│   ├── summary files
│   └── validation outputs
```

---

## Method Summary

Let $\mathcal{A}$ denote the primary archive maintained by RS-CMSA-ESII, and let $\mathcal{C}_f$ denote the feasible and finite passive secondary candidate archive. The final candidate pool is

$$
\mathcal{P}=\mathcal{A}\cup\mathcal{C}_f .
$$

The final reported solution set is obtained using density-filtered score-aware reporting:

$$\mathcal{R}_{\mathrm{DF\mbox{-}SCA}}=\mathcal{D}_{\Delta_f,\tau_\rho}\left(\mathcal{A}\cup\mathcal{C}_f\right),$$

where $\mathcal{D}_{\Delta_f,{\tau_\rho}}(\cdot)$ denotes objective-value filtering and density-based duplicate control.

For a secondary candidate $x$, its normalized distance to the current reported set is

$$\rho(x)=\min_{y\in\mathcal{R}}\left\|\frac{x-y}{u-l}\right\|_2 ,$$

where $u$ and $l$ are the upper and lower bound vectors. The candidate is inserted if

$$
\rho(x)\geq \tau_\rho(D).
$$

The density threshold is dimension-scaled:

$$
\tau_\rho(D)=\alpha_\rho\sqrt{D}.
$$

The final selected setting is

$$
\Delta_f=2\times 10^{-3},
\qquad
\alpha_\rho=2\times 10^{-3}.
$$

Thus,

$$
\tau_\rho(D)=2\times 10^{-3}\sqrt{D}.
$$

---

## Parameter Settings

| Parameter | Value | Description |
|---|---:|---|
| Base optimizer | RS-CMSA-ESII | Original covariance matrix self-adaptation ES with repelling subpopulations |
| Evaluation budget | $FE_{\max}=20000D$ | Official CEC 2026 setting |
| Secondary archive | Passive | Stores restart-best candidates only |
| Value-window threshold | $\Delta_f=2\times10^{-3}$ | Removes weak secondary candidates |
| Density coefficient | $\alpha_\rho=2\times10^{-3}$ | Final selected density-filtering coefficient |
| Density threshold | $\tau_\rho(D)=\alpha_\rho\sqrt{D}$ | Dimension-scaled duplicate-control threshold |
| Final rule | DF-SCA | Density-filtered score-aware reporting |
| Use of true minima during optimization | No | True minima are used only for offline scoring and analysis |

---

## Tested Development Variants

Several variants were evaluated during development.

| Variant | Description | Final decision |
|---|---|---|
| Primary-Archive Baseline | Original RS-CMSA-ESII primary-archive reporting | Baseline |
| SA-Strict | Secondary archive with strict cleanup | Useful but conservative |
| SA-Relaxed | Secondary archive with relaxed reporting | Higher RPR but lower precision |
| SCA-Medium | Medium score-aware reporting | Strong safe backup |
| DimSCA | Dimension-aware reporting | Positive but less stable |
| NBNC-SCA | Nearest-better-neighbor-style cleanup | Too aggressive |
| CapSCA | Candidate-count cap | Precision improved but RPR dropped |
| DF-SCA | Density-filtered score-aware reporting | Final selected rule |
| AAR-SCA | Archive-aware restart initialization | Not retained |
| TLLS-SCA | Lightweight two-level local refinement | Negligible gain |
| CMAR-SCA | Covariance-shaped local refinement | Negligible gain |

The final selected method is **S-CARD-CMSA with DF-SCA reporting**.

---

## Density-Filtering Sensitivity

The density-filtering coefficient was tested using three settings:

| Variant | $\alpha_\rho$ | Mean RPR | Mean precision | Mean F1 | Mean score | Mean $N_{\mathrm{sol}}$ |
|---|---:|---:|---:|---:|---:|---:|
| DF-SCA-Low | 0.001 | 0.558904 | 0.849994 | 0.643001 | 0.600953 | 10.153 |
| **DF-SCA-Balanced** | **0.002** | **0.558904** | **0.870491** | **0.650922** | **0.604913** | **9.819** |
| DF-SCA-High | 0.005 | 0.558954 | 0.810685 | 0.637740 | 0.598347 | 10.494 |

The balanced setting, $\alpha_\rho=0.002$, was selected because it gave the best RPR-F1 trade-off.

---

## Development and Validation Results

On the 320-run development subset, the final DF-SCA rule improved the score compared with SCA-Medium:

| Rule | Mean RPR | Mean precision | Mean F1 | Mean score | Mean $N_{\mathrm{sol}}$ |
|---|---:|---:|---:|---:|---:|
| SCA-Medium | 0.5589 | 0.8516 | 0.6434 | 0.6012 | 10.14 |
| **DF-SCA** | **0.5589** | **0.8705** | **0.6509** | **0.6049** | **9.82** |

On the broader 768-run validation subset:

| Rule | Mean RPR | Mean precision | Mean F1 | Mean score | Mean $N_{\mathrm{sol}}$ |
|---|---:|---:|---:|---:|---:|
| SCA-Medium | 0.805386 | 0.907076 | 0.827346 | 0.816366 | 13.798 |
| **DF-SCA** | **0.805386** | **0.914775** | **0.831271** | **0.818328** | **13.589** |

The final submission-scale run produced 960 output files corresponding to

$$
16\times 15\times 4 = 960
$$

PID-PIN-dimension combinations.

---

## How to Run

### 1. Prepare the CEC 2026 benchmark files

Place the official CEC 2026 benchmark files in the MATLAB path. The code expects the official `ProblemMM` interface and related benchmark files.

Example:

```matlab
addpath(pwd);
addpath(genpath('path_to_CEC2026_benchmark'));
```

### 2. Run a small test

Edit the main batch file to use a small subset, for example:

```matlab
pidList      = [1];
instanceList = [1];
dimList      = [2];
runList      = 1:1;
```

Then run the corresponding main script, for example:

```matlab
main
```

### 3. Run the full official set

For the final submission-scale output, use:

```matlab
pidList      = 1:16;
instanceList = 1:15;
dimList      = [2 5 10 20];
```

The output should contain 960 CSV files.

---

## Submission File Format

Each submitted CSV file should correspond to one PID-PIN-dimension case and contain:

- Columns $1,\ldots,D$: decision variables;
- Column $D+1$: objective value.

A typical file naming convention is:

```text
pid**pin**dim**.csv
```

Please follow the official CEC 2026 technical report instructions for the exact naming and packaging requirements.

---

## Important Notes

- The algorithm does **not** use true global minimizers during optimization.
- True minima are used only for offline development analysis, scoring, and visualization.
- The passive secondary archive does not influence sampling, covariance adaptation, taboo-region updates, restart logic, or termination criteria.
- SHAP analysis and landscape plots are post-hoc interpretation tools only and are not part of the optimizer.
- The final selected version is the density-filtered setting with $\alpha_\rho=0.002$.

---

## Participant

**Dikshit Chauhan**  
Department of Electrical and Computer Engineering  
National University of Singapore  
Email: dikshitchauhan608@gmail.com

---

## Citation

If you use this code or method, please cite the corresponding method report or repository:

```bibtex
@misc{chauhan2026scardcmsa,
  author       = {Dikshit Chauhan},
  title        = {S-CARD-CMSA: A Score-Aware Candidate Archive with Density-Filtered Reporting for Multimodal Optimization},
  year         = {2026},
  note         = {Method report for the IEEE CEC 2026 Competition on Benchmarking Niching Methods for Multimodal Optimization}
}
```

Also cite the original RS-CMSA-ES and RS-CMSA-ESII papers, since S-CARD-CMSA is built on RS-CMSA-ESII.

---

## References

1. Ahrari, A., Deb, K., and Preuss, M.  
   “Multimodal optimization by covariance matrix self-adaptation evolution strategy with repelling subpopulations.”  
   *Evolutionary Computation*, 25(3), 439--471, 2017.

2. Ahrari, A., Elsayed, S., Sarker, R., Essam, D., and Coello Coello, C. A.  
   “Static and dynamic multimodal optimization by improved covariance matrix self-adaptation evolution strategy with repelling subpopulations.”  
   *IEEE Transactions on Evolutionary Computation*, 26(3), 527--541, 2022.

3. Ahrari, A., Fieldsend, J. E., Preuss, M., Li, X., and Epitropakis, M. G.  
   “Experimental setup for CEC 2026 competition on benchmarking niching methods for multimodal optimization.”  
   Technical Report CIS-TF-MMO-TR-2026-001, IEEE CIS Task Force on Multimodal Optimization, 2026.

---

## License
- MIT License for permissive open-source release.
