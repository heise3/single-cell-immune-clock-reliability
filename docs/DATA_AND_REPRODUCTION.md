# Data and reproduction notes

Seed: 20261025. Data-source accessions and model references are identifiers for retrieval, not a transfer of data permissions.

| Resource | Retrieval source | Role |
|---|---|---|
| Allen Sound Life / HIHA | [GEO GSE271896](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271896); CELLxGENE mirror used after mapping | Main longitudinal cohort |
| Terekhova atlas | [HCA project 896f377c-8e88-463e-82b0-b2a5409d6fe4](https://explore.data.humancellatlas.org/projects/896f377c-8e88-463e-82b0-b2a5409d6fe4) | External longitudinal cohort |
| OneK1K | [Yazar et al., Science](https://doi.org/10.1126/science.abf3041); CELLxGENE Census | Masking validation and model training |
| AIDA | [Kock et al., Cell](https://doi.org/10.1016/j.cell.2025.02.017); CELLxGENE Census | Technical repeats and age validity |
| BCG source data | [Li et al., Nature Aging](https://doi.org/10.1038/s43587-025-00819-z), Fig. 5 source data | 38 paired CD8 measurements |
| Published clock | [scImmuAging](https://github.com/CiiM-Bioinformatics-group/scImmuAging), commit `112f8fab10c72aec6e7a7d83762a40117a1e6e24` | Published coefficients and native validation |

Reproduction requires the pinned data release and annotation, full library-size denominators, gene mapping and frozen selection files. Live retrieval can return different versions; do not silently treat those as the original inputs. The archived analysis recorded input SHA256 values at each stage. This public snapshot contains code only; all empirical source tables, individual-level ledgers and trained models must be restored separately under the relevant terms before data-dependent scripts can run.

The inherited K0 archive was built before this clock analysis. Its metadata-only window code is in `cohort_audit/work/`. To run `outputs/round2/step0_population.py`, set `HIHA_K0_SOURCE` to the restored `HIHA_gate_audit` archive; the default is `cohort_audit/outputs/HIHA_gate_audit`. The script checks 68 donors and 320 eligible visits. This is not permission to reconstruct a different cohort after viewing outcomes.

Main final scripts use mean predictions over 100 pseudocells of 15 cells; the independent A1 masking gate uses each donor's median pseudocell prediction, as specified for that gate. Do not conflate these summaries. Leave-one-donor-out calibration excludes every visit and every bootstrap copy of the held-out original donor. Main bootstrap intervals use 2,000 donor resamples.

For external replication, restore the Round 6 inputs, then use the source in this order: `check_Ka_access.py`, `inspect_B_input_format.R`, `complete_Ka.py`, `check_Kb.py`, `fetch_Kc_metadata.py`, `check_Kc.py`, `download_primary_matrices.py`, `run_count_QA.py`, `run_predictions.py`, `analyze_reliability.py`. Original K failures stop downstream work. The study used 101 repeated donors and nominal one/two-year waves; exact visit dates were not in the retrieved cell metadata. These are not asserted to be exact 365/730-day intervals.

Allen annual pairs span recorded vaccination. Only the matching 75–105-day comparison supports the no-recorded-intervening-vaccination description. Batch and visit variance components are model-conditional and do not identify an exclusively technical cause. The BCG simulation concerns responder selection conditional on the observed cohort mean, not a group-average vaccine effect.

Some old receipt scripts mention a separate proposed P2 preregistration label. That administrative label is not a public preregistration of this manuscript. Local staged, time-stamped plans must not be described as prospective registration on an independent public platform.
