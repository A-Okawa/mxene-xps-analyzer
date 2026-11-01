> **Historical report (v1.2.0).** Current test counts and environment: `release/release_metadata_v1.6.0.json`.

# Test report — v1.2.0 (2026-10-16)

## Environments
| Env | Python | Streamlit | numpy | pandas | scipy | Result |
|---|---|---|---|---|---|---|
| `.venv` from `requirements.txt` | 3.12.2 | 1.66.0 | 2.5.4 | 2.3.3 | 1.18.1 | **97 passed** (with private workbook) / 96 passed + 1 skipped (without) |
| anaconda base | 3.12.2 | 1.32.0 | 1.26.4 | 2.2.3 | 1.13.1 | 96 passed + 1 skipped |

## Handoff §29 acceptance criteria

| § | Criterion | Test(s) | Result |
|---|---|---|---|
| 29.1 | (Nb,Ti)2CTx, LiF/HCl, Al, defects → Survey, Ti 2p, Nb 3d, C 1s, O 1s, F 1s, Cl 2p, Al 2p, Li 1s. Candidates 4/3/4/4/2. Vacancy-peak warning | `test_regions_29_1`, `test_fit_candidate_counts_and_warnings_29_1` | PASS (exact region set) |
| 29.2 | C=1, Ti=0.4, Nb=1.2 → M_norm 1.6, Δ 0.4, Ti_norm 0.4, Nb_norm 1.2 | `test_mixed_metal_deficiency_29_2` | PASS (abs 1e-12) |
| 29.3 | V lattice 0.10 + oxide 0.05 → M unchanged, V_dopant_total_norm 0.15. Enters M only via reclassification | `test_dopant_separation_29_3` | PASS |
| 29.4 | Synthetic Ti 2p: convergence, BE/FWHM recovery, grouped stability, no negative areas, identifiability handling | `test_ti2p_refit_29_4[linear/shirley]`, `test_ti2p_matches_prototype_numbers` | PASS |
| 29.5 | Synthetic Nb 3d refit | `test_nb3d_refit_29_5[linear/shirley]` | PASS |
| 29.6 | small perturbation → PASS. +1 eV → FAIL. missing → CHECK-MISSING. grouped PASS + internal Ti differs → warning | `tests/test_validation.py` (8 tests) | PASS |
| 29.7 | missing framework C, missing RSF, descending BE, header/no header, NaN/non-numeric rows, missing provider components, empty optional components | `test_quantification.py`, `test_io.py`, `test_fitting_synthetic.py`, `test_app.py` | PASS |
| 29.8 / 10.4 | Legacy workbook regression ≤ 1e-10 | `test_legacy_workbook_regression` | **PASS: 186 values (31 columns × 6 samples, incl. QC strings), max abs error = 0.0** |

## User decisions D1–D3 (`tests/test_spec_decisions.py`, 14 tests, all PASS)
* D1: ti3c2_brette Ti 2p equals §27 (and not the §12 table). Brette SI 2p3/2, 2p1/2 and C 1s values are stored as
  literature_reference. Window edits apply, are logged, are validated (start inside window; role/ΔSO not editable),
  leave the KB unchanged, and are honoured by the fit.
* D2: nb2c_halim Nb 3d = 204.8/207.3/209.2 eV, FWHM 1.1/1.7/1.0. Windows/ΔSO/ratio/role unchanged; (Nb,Ti) profile
  unchanged. Not hard-fixed: on a spectrum with Nb_MXene at 204.95 eV the fit returns ≈204.95, not 204.8. Request text
  shows the reference values.
* D3: Halim2016 kept; Chen2024 (year 2024, copyright/online 2023). No stale IDs in KB keys, component sources or code.
* All three decisions are recorded as resolved in `spec_discrepancies`, and the KB version is 1.1.0.

## Brette SI O 1s / F 1s references (`tests/test_brette_o_f_references.py`, 26 tests, all PASS)
* Exact HF10 (S5-B, S8, S4-B) and HF48 (S5-A) centers/FWHMs/SI labels per component, including both HF48 P4 rows.
  No sputtered values are stored, and F_residue has no reference value.
* §27 O 1s / F 1s starting values, windows and roles are unchanged.
* Metadata only: fits are bit-identical with references removed or replaced by absurd values. Validation tables are
  identical. Engine modules never read the field. The field is not editable. Stripping all references from the KB
  leaves every fit model identical. A reference outside the window does not invalidate a component.
* ~531 eV: O_disordered stays `termination_O_disordered`, carries the Brette disorder/T-group-vacancy
  interpretation ("NOT a uniquely assigned OH state"), holds only P2 references, and is quantified separately from OH_like.

## Synthetic refit numbers (fixtures in `tests/fixtures`, nbti profile)

| Spectrum | Background | R² | max abs ΔBE (eV) | max abs ΔFWHM (eV) | framework-role area error | max individual component area error |
|---|---|---|---|---|---|---|
| Ti 2p | linear | 0.99984 | 0.110 | 0.217 | 7.3 % | 34 % |
| Ti 2p | Shirley | 0.99983 | 0.140 | 0.222 | 7.1 % | 38 % |
| Nb 3d | linear | 0.99992 | 0.019 | 0.025 | 2.9 % | 7.7 % |
| Nb 3d | Shirley | 0.99992 | 0.020 | 0.062 | 2.2 % | 16.6 % |

The Ti 2p linear case reproduces the v0.1 prototype numbers (R² 0.99984, 0.110 eV, 0.217 eV, 7.3 %, ~34 %).
The new identifiability diagnostic flags Ti_midBE/Ti_highBE (r = −0.89) and Ti_highBE/Ti_oxide (r = −0.85) for Ti 2p,
and Nb_MXene/Nb2O5 and Nb2O5/NbF5_like (r = −0.93) for Nb 3d. This matches the §19 finding: excellent fit, non-unique
internal decomposition.

Fits on the descending, header-less, whitespace-delimited file with junk rows are identical to the clean CSV
(`test_descending_headerless_input_gives_identical_fit`).

## UI tests
* `tests/test_app.py` (streamlit AppTest, headless): all 5 tabs render, request content present, empty-table
  calculation gives a clear error, Nb2CTx / Ti3C2Tx / Custom render, demo fit runs, and switching profile after a fit
  does not crash (regression for a stale-state bug found during E2E).
* `tools/e2e_screenshots.py` (real Chrome via Playwright): uploads the synthetic peak table → corrected_area →
  Calculate. Uploads a descending/no-header Ti 2p file + provider table → Run fit → grouped PASS / component REVIEW.
  Then KB tab, report tab, and the Ti3C2Tx legacy model with synthetic areas. Screenshots: `docs/screenshots/1…6_*.png`.

## Bugs found and fixed during testing
1. An all-blank peak table raised a confusing internal error → now "no component areas entered".
2. Fit results were keyed by core level only. A Ti 2p fit from the (Nb,Ti) profile was re-plotted against the
   Ti3C2Tx model after a profile switch (KeyError) → keyed by profile + core. Quantification results are also tagged by profile.
3. The §28 template mismatch (contaminant C rows) → aligned with the §28 example.
