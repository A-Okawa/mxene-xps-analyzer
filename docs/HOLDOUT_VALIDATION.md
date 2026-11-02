# Prospective hold-out validation protocol (v1.6.0-rc1, frozen baseline f466b94)

**Status (2026-11-02): prospective hold-out validation — not yet performed.** No independent dataset exists yet.
The only real XPS dataset, `261030_sasakiXPS` (samples 1–5), is the **calibration/development dataset**: it was used to
develop and calibrate the model and is therefore **not** independent validation. The frozen first-pass harness is kept
for future prospective validation and was run on the calibration set only to provide the reference side of a future
comparison.

## Future prospective-validation workflow
1. Obtain a genuinely independent XPS dataset measured after this frozen release candidate.
2. Run `tools/validation/holdout_first_pass.py` against the frozen baseline (`f466b94`: app 1.6.0, KB 1.5.0).
3. Preserve the untouched first-pass outputs (never overwritten; later runs use a new label).
4. Review the A–D classifications before changing the KB.
5. Only after prospective validation decide whether a scientific-model revision is necessary.

## Frozen state (must not change during the first pass)
App 1.6.0 · KB 1.5.0 (`data/knowledge_base.json` sha256 `b5d49f7e…b60`) · baseline commit `f466b94` · all peak models,
constraints, windows, splittings, warning thresholds, quantification rules, assignments and open issues as in that
commit. `tools/validation/holdout_first_pass.py` refuses to run if the app/KB version differs, if the KB file differs
from `f466b94`, or if `src/`/`data/` have uncommitted changes.

## Procedure
1. Put the hold-out provider report (MultiPak .pptx or Curve-Fit-Summary text) and raw exports (`*_multi.csv`, optional
   `*curve.csv`) in one folder; optionally a metadata CSV (template `templates/next_dataset_example.csv`).
2. First pass (existing importers and fitting only; no software fix before recording):
   `python -I tools/validation/holdout_first_pass.py --data HOLDOUT_DIR --label holdout --out private/validation/f466b94 --metadata META.csv`
3. Compare with the calibration set:
   `python -I tools/validation/compare_sets.py private/validation/f466b94/first_pass_calibration.json private/validation/f466b94/first_pass_holdout.json private/validation/f466b94`
4. Review every auto-suggested issue category (A software/import, B numerical, C expected limitation of the frozen
   model, D potentially new chemistry / unsupported component). Only A may be fixed before the full review, and only
   without changing the KB or numerical model; the first-pass JSON is never overwritten (later runs use a new label).

Outputs (private, git-ignored): `private/validation/f466b94/first_pass_<label>.json` (all per-region results,
messages with severity, issues, input sha256, frozen-state record), `issues_<label>.csv`, `comparison.{json,md}`.

## Auto-classification rules (suggestions only)
- A: report/raw import failure, sample assembly failure, structure not inferable, refused fit (wrong core level / span).
- B: optimizer failure, non-finite result, covariance not estimable, sensitivity run failure.
- C: statistical non-identifiability / boundary hits / model-sensitive quantities, unsupported regions, provider
  single-line reports, provider components matching a documented calibration observation.
- D: provider components outside every KB window that are not documented calibration observations.
Patterns are compared across sets by sample-independent signatures (reproduced / new / not seen).

## Calibration/development-set reference (samples 1–5, frozen first pass; not a validation result)
46 fits (raw full region + provider curve window), 100 % fit success, 0 import failures; RMSE median 65.6
(raw full region 179.7, provider window 50.7); provider-vs-Python QC FAIL 15 / CHECK-MISSING 12 / REVIEW 12 /
PASS 7; identifiability non-identifiable 105, weakly 28, identifiable 14, boundary-limited 9; 31 boundary hits;
7 unresolved provider components; 33 provisional quantities; 170 model-sensitive spans; messages WARNING 332 /
CAUTION 101 / INFO 22; auto categories C 327, D 6 (the documented provider-only components O 1s ~529.3 eV and
Al 2p ~71.3 / 75.3–76.1 eV). Software check of the full pipeline: dry run on the synthetic example passes (not a
hold-out result).
