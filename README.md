# Single-cell immune clock repeat-measurement reliability

Analysis code accompanying **Repeat-measurement error limits individual-level inference from single-cell immune aging clocks**.

Authors: Peng Liang, Jie Wu, BaoLong Li, XiangYan Wu and Xiaona Zhu. Peng Liang and Jie Wu share first authorship. Manuscript status: submission preparation; no acceptance or journal endorsement is implied.

This is a source-code snapshot of the staged analysis, without empirical datasets or result tables. It includes the published sc-ImmuAging audit, Allen RNA age metric (RAM) scores, comparison models, batch and sampling variance analyses, BCG responder sensitivity analysis, external longitudinal replication and figure scripts. The fixed seed is **20261025**. Publication packaging did not rerun or change empirical analyses.

## Quick verification without participant data

```bash
python3 -m venv work/analysis_env
work/analysis_env/bin/python -m pip install -r requirements.txt
work/analysis_env/bin/python scripts/verify_release.py
work/analysis_env/bin/python -m unittest discover -s tests -v
```

The requirements are the recorded analysis environment, not a claim that all versions are available on every operating system. Python 3.11 was used for the empirical analyses. R 4.6.1 and glmnet 5.0 were used for native checks; additional packages and exact recorded versions are listed in `environment/`. Several historical controllers assume the virtual-environment path above. Run commands from the repository root. POSIX/macOS or Linux is needed for controllers using `fcntl`.

Synthetic checks cover absolute-agreement ICC, SEM/MDC and donor-excluded calibration, including unequal visit counts and bootstrap copies. They do not establish empirical reproducibility. No full clean-environment reanalysis was run for this publication snapshot.

## Analysis map

| Stage | Main source files | Scope |
|---|---|---|
| Frozen cohort windows | `cohort_audit/work/audit_windows_v2.py`; `outputs/round2/step0_population.py` | Allen K0: 68 donors, 320 visits; retrieve original metadata/archive first |
| Gene mapping and masking | `outputs/round2_A1/map_clock_genes.py`, `run_G1_prime.py` | Independent donor-level masking; A1 supersedes the small official-example test |
| Published clocks and RAM | `outputs/round2_A1/score_main_clocks.py`, `compute_reliability.py` | Predictions, equal-donor-weight LODO, bootstrap ICC and variance components |
| Repair comparisons | `outputs/round3/prepare_training_design.py`, `train_repair_models.py`, `evaluate_repair_models.py` | OneK1K training; variants a–e and holdout accuracy |
| Full-cohort comparisons | `outputs/round4/evaluate_F2.py`, `analyze_technical_reliability.py` | Paired-bootstrap differences, batch analysis and repeated cross-fitting |
| Model/validity comparisons | `outputs/round5/analyze_P1_P2.py`, `analyze_P3.py`, `analyze_P4_technical.py` | Published vs matched models, covariates, variance and AIDA technical repeats |
| External longitudinal replication | `outputs/round6/` | Published sc-ImmuAging only; Terekhova earliest/latest healthy visit pairs |
| Figures | `outputs/manuscript/scripts/` | Original figures and source derivation |
| Current Figure 4 and S1 | `outputs/manuscript/aging_cell_submission/scripts/` | Journal wording and S1 legend at upper right |

Only code, documentation, environment records and synthetic validation are distributed. The numbered directories preserve dependencies from the staged analysis; they are not independent executable examples. Historical search/report/packaging scripts are retained for inspection and can expect local records not distributed here. `train_repair_models_coordinate_descent_historical.py` is historical; `train_repair_models.py` is the final training implementation. Round 1 and the small official-example G1 results are superseded and must not be substituted for the final A1 results.

## Inputs and scope

Raw matrices, participant-level metadata, donor/visit ledgers, individual predictions, trained model artifacts, credentials and manuscript files are **not distributed** in this code repository. Obtain inputs under their original providers' terms; see [data and reproduction notes](docs/DATA_AND_REPRODUCTION.md). Exact empirical reproduction also requires the study's frozen selections, input hashes and trained comparison-model artifacts. This repository is not a self-contained empirical rerun bundle. No empirical CSV tables are included. All figure scripts require the relevant original study inputs or source tables, except the design schematic (Figure 1).

Original upstream model repositories are referenced by commit in `outputs/prepare_models.py`; their code and weights are not vendored here. The original sc-ImmuAging, PHARE and other upstream licenses remain applicable. No new license for this project's original code has yet been selected; public visibility does not by itself grant a general reuse license.

`provenance/source_manifest.json` records original and release SHA256 values. Only the external K0 archive location and optional font lookup were changed for portability; analysis logic was preserved. `SHA256SUMS.txt` hashes the released files. See `CITATION.cff` for software citation metadata. This GitHub URL is a code repository, not a DOI or a public preregistration record.
