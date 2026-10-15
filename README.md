# MAX/MXene XPS Analyzer (v1.2.0)

**Current status: literature-based model frozen; waiting for real XPS calibration.** (tag `v1.2.0-precalibration`, KB v1.2.0 — restart guide: [`docs/RESUME_AFTER_XPS.md`](docs/RESUME_AFTER_XPS.md))

Streamlit lab app for **XPS measurement planning, provider fitting requests, defect/dopant-aware
quantification, and an independent Python refit / QC** of MAX/MXene spectra. The first target is
defect-rich **(Nb,Ti)2CTx**. Also included: Ti3C2Tx (Brette2024 model plus the lab's legacy Excel model) and Nb2CTx.

Specification (source of truth): `CLAUDE_CODE_XPS_APP_HANDOFF.md`. Prototype v0.1 was used only as a behavioural reference.

> The app does **not** map binding energy to a unique chemical state. Composition change, M-site
> deficiency, termination disorder, dopant incorporation, oxidation and secondary phases are kept as separate
> variables. XPS-only results are labelled *apparent* / *candidate*.

## Workflow (one tab per step)

| Tab | What it does |
|---|---|
| Sidebar | Sample setup (structure, M/X targets, M elements, etchants, A element, defects, dopants, nominal M stoichiometry) and measurement metadata (X-ray source, pass energy, step, BE reference, sputtering) |
| 1. Request generator | Region list (§7 rules) plus an editable, downloadable fitting request containing every candidate component's start BE, window, FWHM, ΔSO, area ratio, role and source |
| 2. Peak table & quantification | Template, upload, and editable table with a role dropdown. Role changes are logged as user reclassifications. Supports raw/corrected/atomic-% input; RSFs are required in raw mode (no silent 1.0). Carbide-C normalization, apparent M deficiency, per-metal nominal-vs-observed fractions, O/F/Cl metrics (raw and corrected), dopants by role, orthogonal-evidence record. For Ti3C2Tx you can choose the **legacy workbook model** |
| 3. Provider ↔ Python validation | Raw spectrum import (any order, header optional), Shirley/linear background, pseudo-Voigt doublet fit with NNLS start, R²/RMSE/AIC/BIC, residual plot (BE axis inverted), identifiability diagnostics, component-level and role-grouped QC with PASS/REVIEW/FAIL/CHECK/CHECK-MISSING |
| 4. Knowledge base | Principles, profile notes, every component with its caveat, source key and DOI, role ontology, the Han solid-solution table, open spec discrepancies, JSON download |
| 5. Report & export | Project/session JSON (setup, metadata, KB version, RSFs, tables, results, warnings, fit settings, QC thresholds, sources) and a consolidated HTML report |

Screenshots: `docs/screenshots/` (synthetic data only).

## Run locally (macOS / Linux / Windows)

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

## Deploy for students

**Data rule:** unpublished raw spectra and the historical Excel workbook must never be committed. Data are
uploaded at runtime only. `.gitignore` blocks `private/`, `*.xlsx`, and non-synthetic files. No app feature
sends data to an external AI/API service.

* **Streamlit Community Cloud:** push this repository (code + `data/knowledge_base.json` + synthetic
  fixtures only) to GitHub → *New app* → main file `app.py`. Uploaded spectra then go to Streamlit's servers,
  so for unpublished data use a private repo/app with viewer restrictions, or one of the options below.
* **Institutional / private server:** `streamlit run app.py --server.port 8501 --server.address 0.0.0.0`
  behind the lab VPN or a reverse proxy with authentication.
* **Local:** students run the commands above. Nothing leaves the PC.

## Tests

```bash
pip install pytest
python -m pytest -q                                   # 96 pass, legacy regression skipped
XPS_LEGACY_WORKBOOK=/path/to/mxene_surface_functional_group_quant_template_v2.xlsx python -m pytest -q   # 97 pass
```

The legacy regression reads the private workbook's cached Excel results and checks all 31 output columns × 6
samples against the Python legacy model (tolerance 1e-10). You can also copy the workbook to
`private/legacy_workbook.xlsx` (git-ignored).

Browser E2E + screenshots (dev only; needs `pip install playwright` and Google Chrome):
```bash
streamlit run app.py --server.headless true --server.port 8599 &
python tools/e2e_screenshots.py
```

## Repository layout

```text
app.py                      Streamlit UI only (no scientific logic)
src/mxene_xps/
  models.py                 role ontology, status vocabulary, typed schemas (FitComponent validates every KB entry)
  knowledge_base.py         versioned KB loader + validation (every component needs source + caveat)
  request_generator.py      §7 region rules, §8 request text, §28 peak template
  quantification.py         §9 role-based engine, §10 exact legacy workbook model, reclassification audit
  fitting.py                §17 Shirley/linear, pseudo-Voigt doublets, NNLS init, metrics, identifiability diagnostics
  validation.py             §18 component + role-grouped QC and summary logic
  reporting.py              §21/§22 context warnings, evidence summary, project JSON, HTML report
  io.py                     raw-spectrum importer registry, peak-table validation
data/knowledge_base.json    KB v1.0.0 (built by tools/build_kb.py from the handoff §27 seed)
tools/                      build_kb.py, kb_seed_from_handoff_s27.json, make_fixtures.py, e2e_screenshots.py
tests/                      pytest suite (§29) + synthetic fixtures
docs/                       scientific_basis.md, references.md, TEST_REPORT.md, screenshots/
```

## Peak-table schema

Required: `sample, core_level, component, area, role`. Recommended: `element, BE, FWHM, notes`.
A blank area means "not fitted" and is listed in the warnings. A non-numeric area, a negative area or an unknown role
is rejected. Roles: see the *Role ontology* table in tab 4 or `src/mxene_xps/models.py`.
