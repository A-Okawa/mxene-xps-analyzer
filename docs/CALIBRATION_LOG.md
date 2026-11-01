# Calibration log (branch `calibration/real-xps`)

Starting point: tag `v1.2.0-precalibration`. First real data: provider (PHI MultiPak) fits of Nb2AlC,
(Nb0.9Ti0.1)2AlC, (Nb0.6Ti0.4)2AlC, Nb2CTx and (Nb0.6Ti0.4)2CTx (private; never committed).

## Phase 1 — I/O only (app 1.2.1, KB unchanged at 1.2.0)

| Problem found on real data | Fix |
|---|---|
| `*_multi.csv` (9 regions) was read as ONE spectrum (−4.6…700 eV); an Nb 3d fit "succeeded" with R² 0.94 | `parse_two_column` refuses labelled multi-block files with BE discontinuities; new `parse_phi_multiregion` splits them into one spectrum per region (UI: region picker) |
| No guard against implausibly wide ranges | `check_single_region`: hard stop above 60 eV span, warning above 40 eV / duplicate BE; enforced inside `fit_core_level` (`max_span_eV`) |
| Provider table only available as MultiPak "Curve Fit Summary" text in a .pptx; spin-orbit lines reported as separate bands | `mxene_xps.multipak`: parses summaries (.pptx or text), pairs lines by region, splitting, area ratio and proximity (importer heuristics, not KB values); lone lines are flagged, totals estimated and marked |
| Provider labels were lost when renamed to KB names | Three separate layers: raw bands (untouched) → paired rows (carry provider labels) → tentative KB mapping (method `auto_window_tentative` / `user`, logged). Unmapped rows keep the provider label, role `other` |
| Fit window could not exclude neighbouring lines | Optional user fit window (`crop`), recorded in fit settings |

Tests: `tests/test_phi_multipak_io.py` (synthetic; optional real-data check with `XPS_SASAKI_DIR`).
