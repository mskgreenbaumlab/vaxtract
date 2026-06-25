# figure_benchmark — RESULTS (Rojas, PMID 37165196)

- **Date:** 2026-06-03
- **Run:** vision arm, 2 curated panels, `claude-opus-4-8` via `claude -p`.
- **Cost:** Fig 1G $0.45, Fig 4B $0.28 — **~$0.72 total**. Wall ~21s / ~11s.
- **Ground truth:** Nature Source Data xlsx (machine-exact). Scoring is value-based
  (`by_nearest`, single-series collapsed — the vision model invents its own series
  labels, so matching is by value, not by a label it could never reproduce).
- Report artifacts: `outputs/figure_benchmark/report.{json,html}` (gitignored).

## Per-panel metrics (vision vs source-data truth)

| Panel | Chart | matched | missed | extra | median \|log10\| | p90 \|log10\| | peak log-RMSE | mean log-RMSE |
|-------|-------|--------:|-------:|------:|----------------:|-------------:|--------------:|--------------:|
| **Fig 4B** | 3 linear bars (% blood T cells) | 3/3 | 0 | 0 | **0.005** | 0.019 | 0.005 | 0.006 |
| **Fig 1G** | dense IFNγ ELISpot scatter, log/segmented y, ~29 dots | 25/29 | 4 | 0 | **0.258** | 0.621 | **0.024** | 0.104 |

`|log10|` is the absolute log10 ratio of read-vs-truth value. 0.301 = 2×, 0.097 ≈ 1.25×.

## Reading the numbers

- **Bars / simple charts (Fig 4B): vision is essentially exact.** Median per-point
  error ~1% (0.005 log10); all 3 bars recovered. Truth `[0.7, 1.7, 2.5]` → read
  `[0.7, 1.6, 2.5]`.
- **Dense log-scatter (Fig 1G): per-series summary excellent, per-point ~2× scatter.**
  - The **peak / maximum magnitude is read almost exactly** (peak log-RMSE 0.024 ≈ 6%;
    truth max 2459 → read ~2500). This is the summary the DB actually stores.
  - **Individual overlapping points carry ~1.8× median error** (0.258) and up to ~4×
    at p90 (0.621); 4 of 29 dots were missed (expected for overlapping marks on a
    small, segmented-axis scatter).
  - Vision correctly tracked the segmented axis (linear 0–500 then compressed
    500/1.5k/3k) — the high points were not catastrophically mis-read.

## Verdict — go/no-go on building a real `read_figure`

**GO, scoped by chart type and readout:**

1. **Bars, small categorical charts, line charts with few points → use vision values
   directly** (well within any reasonable tolerance; ~1% error here). A real
   `read_figure` can emit these as `value=<read>` with normal provenance.
2. **Dense scatter / log panels → use vision for the per-series *summary* (peak/mean
   magnitude), NOT exact per-point values.** Emit the summary as
   `value=<read>, tier='reported', needs_review=true`, and keep per-point dots as
   `needs_review` estimates. The ~2× per-point scatter means exact per-dot values
   should not be trusted unattended.
3. **Proposed tolerance to confirm with Samuel:** accept a vision read when
   `|log10(read/truth)| <= 0.1` (~1.25×) for simple charts; treat dense-scatter
   *summary* reads as acceptable at `peak log-RMSE <= 0.05` (~12%). Per-point dense
   reads stay review-gated regardless.

## Caveats / scope of this run

- **2 panels, 1 paper.** A real go/no-go for production should add a few more papers
  and chart types (KM curves, heatmaps were out of scope — no clean source-data match
  in Rojas).
- **Vector comparison arm not exercised.** Built + unit-tested, but none of Rojas's
  curated panels suit it (contingency = no axis, bars = not circular marks, Fig 1G =
  segmented axis that the linear/log calibration can't model faithfully). Exercising
  it needs a scatter with a clean linear/log axis.
- **Fig 2A dropped:** the source sheet "Figure 2A bar graph" (2×2 contingency) does
  not match the published Fig 2 panel-a line charts; not locatable at that letter.
- **Next step (separate, gated on this):** wire a numeric `read_figure` into the agent
  per verdict items 1–2, keeping dense per-point reads `needs_review`.
