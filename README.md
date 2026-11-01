# MAX/MXene XPS Analyzer — v1.6.0 (release candidate)

App **1.6.0**, knowledge base (KB) **1.5.0**, branch `calibration/real-xps`.
Calibration history: [`docs/CALIBRATION_LOG.md`](docs/CALIBRATION_LOG.md). Machine-readable release record:
[`release/release_metadata_v1.6.0.json`](release/release_metadata_v1.6.0.json).

A local Streamlit app for XPS work on MAX phases and MXenes: measurement/fitting requests, import of raw spectra and
provider fits, an independent Python refit with QC, defect-aware quantification, and cross-sample comparison.
The app does **not** map binding energy to a unique chemical state. Composition, M-site deficiency, termination
disorder, oxidation, dopants and secondary phases are kept as separate variables, and XPS-only results are labelled
*apparent*, *candidate* or *provisional*.

## Three levels that are not equivalent

| Level | Meaning | Where shown |
|---|---|---|
| **1. Fit succeeded** | The optimizer converged; residuals, R², AIC/BIC are reported. Says nothing about uniqueness or chemistry. | tab 3 metrics; `fit_succeeded` in reports |
| **2. Statistically identifiable** | Within this model, the component's parameters are determined by the data (status: identifiable / weakly identifiable / boundary-limited / non-identifiable). Still not a chemical assignment. | tab 3 identifiability table; reports |
| **3. Chemically validated** | The assignment is confirmed by orthogonal evidence (e.g. XRD, XAS, EELS, reference spectra). **The app never claims this**; reports always state `chemical interpretation validated = False`. | reports |

A converged fit can contain non-identifiable components, and an identifiable component can still be chemically
mis-assigned.

## Supported structures, profiles and core levels

| Structure | Profile | Fitted core levels | Basis of starting values |
|---|---|---|---|
| Ti3C2Tx | `ti3c2_brette` | Ti 2p, C 1s, O 1s, F 1s | literature (Brette 2024); legacy workbook model optional |
| Nb2CTx | `nb2c_halim` | Nb 3d, C 1s, O 1s, F 1s | literature (Halim 2018) + calibration candidates |
| (Nb,Ti)2CTx | `nbti_m2c_defect_aware` | Ti 2p, Nb 3d, C 1s, O 1s, F 1s | literature/transferred + calibration candidates |
| Nb2AlC | `nb2alc_max` | Nb 3d, C 1s, O 1s, Al 2p | literature (Su 2024 SI; Näslund 2020) + transferred |
| (Nb,Ti)2AlC | `nbti_211_max` | Ti 2p, Nb 3d, C 1s, O 1s, Al 2p | **transferred** (no direct (Nb,Ti)2AlC peak table) |

**Recognized but unsupported:** Cl 2p, S 2p, Li 1s. These can be imported, displayed and exported, and their provider
peaks can be imported (spin-orbit pairs are kept), but they have **no validated KB fit model**. Fitting any spectrum
with an unrelated core-level model is refused.

## Launch

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py                                   # http://localhost:8501
```

Nothing is sent to external services. Unpublished data are uploaded at runtime only and are never written to the
repository (`.gitignore` blocks data formats; only synthetic fixtures, templates and the synthetic example are tracked).

## Inputs

| Input | Format | Notes |
|---|---|---|
| Single-region spectrum | two columns (BE, intensity); any delimiter, header optional, any order | refused if it contains several labelled blocks with BE jumps; refused (INVALID) for a single-region fit if the BE span is > 60 eV, warned > 40 eV |
| **PHI multi-region CSV** (`*_multi.csv`, `*_su.csv`) | MultiPak ASCII export: `Full` / blank / `<region>` / `1` / data, repeated | split into one spectrum per region; pick the region in tab 3; surveys are flagged "not for peak fitting" |
| **MultiPak provider fit** | `.pptx` report or plain-text "Curve Fit Summary" (+ optional "Atomic Concentration Table") | raw bands kept as reported; spin-orbit lines paired by region, splitting, area ratio and proximity; lone lines flagged (totals estimated); tentative, editable mapping to KB components (edits logged) |
| Generic provider table | CSV `component, BE, FWHM, area[, role, core_level]` | doublet area convention selectable |
| Peak table for quantification | CSV `sample, core_level, component, area, role` (+ element, BE, FWHM, notes) | raw area (RSF required), corrected area, or atomic % |
| **Sample metadata** | CSV, template in `templates/sample_metadata_template.csv` (example for the next five samples: `templates/next_dataset_example.csv`) | sample ID, composition, MAX/MXene, synthesis, etchant, treatment, temperature, date, BE reference, notes, sputtering, comments; dataset-level only, never part of the KB |

## Workflow (tabs)

1. **Request generator** — region list and a fitting request with every candidate component (start BE, window, FWHM,
   ΔSO and its basis, role, source, literature/transferred references, calibration status).
2. **Peak table & quantification** — carbide-C normalization, M_norm, apparent δ, per-metal fractions, termination
   quantities, and the **provisional table**: nominal (calibration candidate excluded from the framework) vs
   alternative (candidate counted as framework), shown side by side and **never averaged**, with the reason each
   value is provisional.
3. **Provider ↔ Python validation** — independent pseudo-Voigt refit (Shirley/linear background), provider-vs-Python
   component and role-grouped QC, the **identifiability table** (BE/FWHM/area with standard errors, bound proximity,
   strongest correlations, anticorrelation partner, stability across alternatives, status), and **fit-sensitivity
   comparison** over predefined alternatives (background, model-window crop, D4 ΔSO grid for Nb 3d only). ΔAIC/ΔBIC are
   shown only where the data points are identical; area-fraction or δ spans ≥ 0.05 are labelled **model-sensitive**.
   The app never chooses the "best chemistry" by AIC/BIC.
4. **Knowledge base** — every component with its role, caveat, literature / transferred references, provider
   observations (display only), calibration status and DOI.
5. **Report & export** — project JSON, HTML report, and the **calibration summary** (JSON/HTML).
6. **Cross-sample comparison** — provider report + optional raw multi-region files + optional metadata → BE and area
   fraction by component, O 1s / F 1s fractions, survey at.% (Al residual, Cl, S, Li where measured), provisional
   quantities and unresolved components; group by phase, composition or treatment.

## Provenance

Every exported value carries one of: `raw_measurement`, `provider_raw_fit`, `provider_paired_doublet`, `python_fit`,
`literature_reference`, `transferred_reference`, `calibration_observation`, `user_override`, `derived_quantity`.
Provider fit values of real data are stored only as calibration observations; the KB refuses to load if an uncited
source is used as a reference or as the starting value of a non-candidate component.

## Messages

| Severity | Used for |
|---|---|
| INFO | neutral notes (import details, inferred metadata, how a value was computed) |
| CAUTION | result usable but chemistry uncertain or provisional (calibration candidates, provisional δ, unresolved O, transferred models, referencing caveats, estimated areas, user edits) |
| WARNING | statistical/model problem limiting interpretation (strong (anti)correlation, boundary-limited or non-identifiable component, model-sensitive quantity, QC FAIL, BE discontinuities) |
| INVALID | refused or failed (wrong core-level model, implausibly wide BE span, failed import/fit, invalid input) |

## Exports

Project JSON (setup, metadata, release block, tables with provenance, fits incl. provider raw bands / pairs / mapping /
user overrides / identifiability / sensitivity, provisional quantities, warnings with severity), HTML report, calibration
summary (dataset, app + KB version, git commit, Python and package versions, fit models, QC, identifiability,
model-sensitive quantities, unresolved components, overrides, references, open issues), and CSVs (fit components,
curves, QC, quantification, cross-sample long table) with provenance columns.

## Example dataset (synthetic)

`examples/synthetic_demo/` — **synthetic** files for onboarding and software validation (single-region, PHI
multi-region, MultiPak text report with atomic-concentration tables, generic provider table, metadata). Not measurements;
not derived from any real data. See its README.

## Tests

```bash
pip install pytest
python -m pytest -q                                            # runs without any private data (private tests skip)
XPS_SASAKI_DIR=/path/to/private/provider_folder \
XPS_LEGACY_WORKBOOK=/path/to/private/legacy_workbook.xlsx python -m pytest -q   # adds private regression tests
```

Browser E2E (synthetic): start `streamlit run app.py --server.headless true --server.port 8599`, then
`python tools/e2e_screenshots.py OUT_DIR` (Playwright + Chrome).

## Scientific limitations (current; not solved in this release)

- The Nb 3d **205.7–205.9 eV component** (`Nb_suboxide_like`) is **not chemically validated**; the O 1s cross-check was
  inconclusive. It is excluded from the framework by default.
- The **apparent M-site deficiency δ is model-dependent and provisional** (candidate excluded vs counted; provider vs
  Python decompositions differ).
- The **Nb 3d spin-orbit splitting (D4) is calibration-sensitive**; the KB constraint stays 2.8 eV.
- **Ti 2p in mixed Nb/Ti samples is affected by Nb 3s overlap**; interpret Ti 2p3/2 first.
- In the MAX profiles, **Nb_MAX and Nb_MAX_2nd are strongly correlated** and should normally be interpreted by their sum.
- In the first real dataset, **sample 2's provider Ti quantification is biased** because only Ti 2p3/2 was fitted.
- The **sample-1 low-BE Nb feature** (`unresolved_low_BE_Nb`) is real extra intensity but **non-identifiable**; it is not
  assigned to any species.
- **Cl 2p, S 2p and Li 1s have no validated KB fit model.**
- **Transferred models** (e.g. `nbti_211_max`) are **not equivalent to direct literature validation**.
- XPS alone cannot give absolute vacancy concentrations; orthogonal evidence is required.
