# Scientific basis and implementation decisions

Source of truth: `CLAUDE_CODE_XPS_APP_HANDOFF.md` (the handoff). This file maps its rules to the code and lists
every place where the implementation had to make a choice the handoff does not settle.

## 1. How the non-negotiable rules (§2, §31) are enforced

| Rule | Where it is enforced |
|---|---|
| No BE → unique state lookup | No code path assigns chemistry from BE. Roles come only from the KB starting model or from an explicit user edit, which is logged (`role_reclassifications`) and exported |
| No absolute vacancy from XPS | Output is named `apparent_M_deficiency_delta`, with the §30 warning. Orthogonal evidence is recorded and displayed, never scored |
| Never force Nb+Ti = 2 | `quantify_sample` normalizes by framework C only (`scale = X_target / C_framework`). The test `test_metal_sum_not_forced_to_ideal` guards this |
| Fit metals independently | Each core level is fitted separately. No shared fractions exist in the model |
| Dopant ≠ substitution | Dopant roles are summed separately per element and role. They enter M only if a row is reclassified to `framework_metal` |
| Mixed-metal O 1s ambiguous | `O_low` / `O_mid` have roles `O_unresolved_*` and `assignment_status = "unresolved"`. They are never counted as Tx, and no UI switch counts them; only a logged per-row role change can |
| F residue ≠ Tx | `residue_F` is excluded unless the user overrides, which triggers a warning |
| Sputtering / BE reference | Measurement metadata in the sidebar. `metadata_warnings` flags unknown/non-Fermi-edge referencing, adventitious-C forcing, and sputtered data |
| No silent RSF = 1 | Raw-area mode raises if any core level lacks a positive RSF. Legacy RSFs are prefilled only on an explicit opt-in |
| No unique Nb vacancy peak / no Ti²⁺/³⁺/⁴⁺ primary model | Guarded by KB tests. Legacy component names are shown with a "historical labels, not oxidation states" notice |

## 2. Knowledge base (v1.2.0)

`data/knowledge_base.json` is built by `tools/build_kb.py` from `tools/kb_seed_from_handoff_s27.json`, which is the
handoff §27 JSON extracted verbatim. `test_numerical_values_identical_to_handoff_seed` checks that every center,
window, FWHM, split, ratio, role and source matches the seed exactly, **except** the six values changed by
decision D2 (listed explicitly in the test).

Each component carries three kinds of numbers, kept apart:

| Field | Meaning | Editable in UI |
|---|---|---|
| `center`, `fwhm` | starting value for the fit | yes (logged) |
| `bounds`, `fwhm_bounds` | implementation fitting windows (software, not literature constants) | yes (logged) |
| `literature_reference` | list of values reported in the cited SI (one entry per SI table/dataset); display/request only, **never a constraint** | no |
| `literature_interpretation` | the cited paper's own interpretation (e.g. Brette ~531 eV O 1s) | no |

User edits to starting values and windows (tab 3) are applied by `knowledge_base.with_overrides`. The KB itself is
not modified. Each edit is recorded as `model_overrides` in the fit settings and the project export. Roles, kind,
ΔSO, area ratio, sources and literature references cannot be edited there.

What was added in addition to the seed (metadata only): `interpretation` caveats where missing, using wording
already in the handoff §12–§16, with no new chemical assignment; `assignment_status`; `provenance`;
`kb_version`/`created_at`/`reviewed_by`/`changelog`; references `Chen2024` and `LegacyExcel`.

## 3. Specification discrepancies — resolved by user decision (2026-10-16)

* **D1 – Ti3C2Tx Ti 2p.** §27 JSON adopted for `ti3c2_brette` (§12 table values not used). Brette SI experimental
  values (HF10 %, before sputtering; handoff §12) are stored as `literature_reference`:
  * Ti 2p3/2: 454.9 (1.0) / 455.6 (1.6) / 456.8 (1.7) / 459.3 (1.8) eV
  * Ti 2p1/2: 460.8 (1.1) / 461.5 (1.6) / 462.4 (1.6) / 464.8 (3.2) eV
  * C 1s: 281.8 (0.6) / 284.1 (1.9) / 286.4 (1.8) / 288.9 (2.0) eV
  
* **Correction (2026-10-16) – Brette O 1s / F 1s.** The Brette SI (`downloadSupplement.pdf`) does give exact
  values. I transcribed them and checked them against the rendered tables. They are stored as `literature_reference`
  entries labelled by SI table, dataset and sample state, and are **reference metadata only**:

  | KB component (Brette label) | HF10 S5-B | HF10 S8 (before sputt.) | HF48 S5-A |
  |---|---|---|---|
  | O_ideal (P1) | 529.7 (1.0) | 529.6 (0.9) | 529.8 (1.0) |
  | O_disordered (P2) | 531.0 (2.0) | 531.1 (2.0) | 531.0 (2.0) |
  | OH_like (P3) | 532.5 (1.6) | 532.5 (1.5) | 532.0 (1.7) |
  | O_adsorbate (P4) | 534.0 (1.6) H2O+OR | 534.0 (2.1) H2O+OR | 533.8 (1.7) OR; 535.8 (2.0) H2O |
  | F_MXene | 684.9 (1.3) | 684.8 (1.3) (Table S4-B) | 685.1 (1.6) |

  Components are mapped using the paper's own P1–P4 labels. After-sputtering rows are not stored (sputtered
  surfaces are not pooled, §22). Brette reports a single F 1s contribution, so `F_residue` has no reference value.
  The §27 starting values and fitting windows are unchanged (e.g. OH_like still starts at 532.2 eV).
  `fitting.py`, `validation.py`, `quantification.py` and `io.py` never read `literature_reference`, and a test enforces this.
* **~531 eV O 1s (Brette interpretation).** `O_disordered` keeps role `termination_O_disordered`, which is
  quantified as `z_O_disordered` separately from `z_OH_like`. It carries `literature_interpretation`: O on disordered
  surfaces (mix and/or vacancies on T-groups), whose magnitude probes vacancy-induced surface disorder; it is **not a
  uniquely assigned OH state**. The OH-like contribution is the separate ~532 eV P3.
* **D2 – Nb2CTx Nb 3d.** `nb2c_halim` starting center/FWHM = Halim SI Nb2CTx values: 204.8/207.3/209.2 eV,
  FWHM 1.1/1.7/1.0 eV (C–Nb–Tx/Nb2O5/NbF5). These are also stored as `literature_reference`. They are
  initial/reference values: the §27 windows (unchanged; they contain the new values) let the fit move.
  The (Nb,Ti)2CTx profile is unchanged.
* **D3 – Source keys.** `Halim2016` kept; Chen standardized to `Chen2024` (Small, year 2024;
  `copyright_or_online_year: 2023`).

Regression tests: `tests/test_spec_decisions.py`, `tests/test_brette_o_f_references.py`.

## 4. Implementation choices where the handoff is silent

* **§28 template:** like the §28 example, `contaminant` C 1s rows are left out of the quantification template (they never
  enter quantification). They remain in the fit model and the request text.
* **Metal-QC tolerance** `|M_target − M_norm| ≤ 0.35` (`metal_QC`, deficiency/excess warnings): this value comes from
  the prototype. It is a software default and user-editable. The legacy model keeps its own ±0.5 rule.
* **Negative Δ (M excess)** gets a neutral warning (possible C under-count / RSF / misclassification). It is not
  interpreted structurally.
* **Validation:** ΔBE is signed (Python − provider, §18.1) and graded on |ΔBE|. A metric the provider did not supply
  (e.g. no FWHM) is graded `CHECK`, not silently PASS as in v0.1. Severity order for the overall row status:
  FAIL > CHECK-MISSING > REVIEW > CHECK > PASS. Components present only in the provider table or only in the Python
  fit are `CHECK-MISSING`.
* **Provider doublet-area convention:** the user states whether provider doublet areas are totals or primary-line
  only. The request text asks the provider to say which.
* **Grouped QC with a provider role** that differs from the KB role uses the provider/user role and reports the
  disagreement.
* **Identifiability diagnostics (§19, §29.4):** the covariance from the least-squares Jacobian gives area–area
  correlations. |r| ≥ 0.8 (editable) is flagged "not uniquely identifiable – prefer role-grouped totals". Parameters
  at a window edge and components below 1 % of the envelope are also flagged. These flags are software heuristics.
* **Shirley:** iterative, trapezoidal cumulative integral (works for non-uniform steps), endpoint averaging optional
  (default 1 point = prototype behaviour). Not calibrated against CasaXPS.
* **Legacy model:** blank inputs behave like Excel (blank = 0). Only framework C must be positive. Operation order
  copies the workbook formulas.

## 5. Assumptions that still need calibration with real historical data (§33)

1. Shirley implementation vs the provider's/CasaXPS background (endpoints, averaging, iteration).
2. Line shape: pseudo-Voigt with free η and shared FWHM/η within a doublet vs the provider's line shapes
   (e.g. GL(30), LA, asymmetric tails for metallic/carbide components).
3. Fixed ΔSO and area ratios (Ti 2p 5.6–5.8 eV / 2:1, Nb 3d 2.8 eV / 1.5:1). Real Ti 2p often needs a wider 2p1/2 FWHM.
4. BE windows and FWHM bounds (software windows) for this provider's instrument and referencing.
5. QC thresholds (0.30/0.50 eV, 25/50 %, 0.10/0.20) and the 0.8 correlation flag.
6. RSFs: instrument/software specific. The legacy RSFs (C 0.314, Ti 2.077, O 0.733, F 1.0, Cl 0.954) are not valid for
   Nb 3d or for another instrument.
7. The 0.6667 F-oxide factor: legacy Ti3C2Tx only. Not transferred to (Nb,Ti).
8. The provider's export format: `io.SPECTRUM_IMPORTERS` currently holds only a generic two-column reader.
9. Mixed-metal C 1s: one `C_MXene` component (282.15 eV, Han2020). Splitting Ti–C vs Nb–C is deliberately not
   implemented.
10. The (Nb,Ti)2CTx profile is a composite (Han + Halim + Brette). No single reference validates it.
