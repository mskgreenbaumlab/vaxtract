# `scale/` — HPC batch extraction (subscription lane)

Run the cancer-vaccine extractor over a corpus of papers on SLURM, one task per PMID.

## The one thing to understand first
This is **not a compute job** — each task is a thin client that mostly sits idle waiting on the
Claude API. So:
- it asks for **tiny CPU/mem + long walltime**, not GPUs/many cores;
- the scaling ceiling is the **subscription's rolling-window quota**, not the cluster;
- therefore **concurrency is pinned to `-j 1`** (one `claude` login can't be shared across parallel
  jobs, and its OAuth token refresh races). This is a serial, resumable overnight queue.

Realistic throughput on a subscription is a handful-to-dozens of papers per usage window, then it
pauses until the window resets and you re-run the same command (it resumes — see below). For true
parallelism you need the paid lane (`auth: api_key` / Bedrock) — built-in, just not funded yet.

## Layout
```
scale/
  Snakefile              one `extract` task per PMID; output = a terminal marker
  config.yaml            EDIT THIS — absolute paths, auth, creds, concurrency
  extract_one.py         per-paper wrapper: run -> validate -> route -> ledger -> marker
  corpus/pmids.txt       the corpus (one PMID per line)
  singularity/           runtime image recipe (Node+claude CLI + Python+deps)
  probe_egress.sh        pre-flight: prove a compute node has egress + a working login
```
Results land under `results:` from config:
```
results/extraction/
  extracted/<PMID>.json     valid records (the deliverable; silver — still needs human sign-off)
  needs_review/<PMID>.json  records that FAILED schema validation (a human looks at these)
  markers/<PMID>.done       terminal marker => never re-extract this paper
  logs/<PMID>.log           full runner transcript per paper
  ledger.csv                scaled RUNS.md: timestamp,pmid,status,subtype,turns,cost_usd,valid,...
```

## Resumability + quota safety (the marker contract)
`extract_one.py` writes the Snakemake rule output (the marker) **only for terminal outcomes**:

| outcome | JSON | marker | exit | next `snakemake` run |
|---|---|---|---|---|
| `ok` | `extracted/` | yes | 0 | skipped (done) |
| `invalid` | `needs_review/` | yes | 0 | skipped (human reviews) |
| `usage_limit` | — | no | 1 | **re-extracted** (after window resets) |
| `error` (net/crash) | — | no | 1 | **re-extracted** |

So when the quota stops you mid-corpus, just re-run the same command after the window resets and it
picks up exactly where it left off — without re-spending quota on finished papers or endlessly
retrying a genuinely-bad one.

## One-time setup
1. **Log in once** on the login node so the CLI stores its subscription creds:
   `claude login` → creds land in `creds_dir` (default `~/.claude`). Ensure that dir is on **shared
   home** so compute nodes can read it.
2. **Build the image** (on a host with internet + fakeroot):
   `singularity build --fakeroot scale/singularity/antvac_extractor.sif scale/singularity/antvac_extractor.def`
   then copy it to the `sif:` path in `config.yaml`.
3. **Edit `scale/config.yaml`** — set `repo`, `sif`, `data_root`, `results`, `creds_dir`.
4. **Stage sources**: each PMID needs `data_root/<PMID>/{main.pdf, supps/*.pdf,*.xlsx}`.

## Pre-flight (do this before the first batch)
```
srun -A your-slurm-account -p componc_cpu -t 00:10:00 --mem 4G \
  bash scale/probe_egress.sh <sif> <creds_dir>
```
Exit 0 ⇒ a compute node has API egress **and** the subscription login works inside the container.

## Run
```
snakemake -s scale/Snakefile --configfile scale/config.yaml \
  --workflow-profile profiles/workflow_profiles/snakemakes/slurmMinimal -j 1
```
Re-run the same line after each usage window until `results/extraction/markers/` has every PMID.
Watch progress with the ledger: `column -s, -t results/extraction/ledger.csv | less -S`.

## Flipping to the paid lane later (no code changes)
In `config.yaml`: set `auth: api_key`, provide `ANTHROPIC_API_KEY` in the job env (or set
`CLAUDE_CODE_USE_BEDROCK=1` + IAM for Bedrock), and raise `max_concurrent` / `-j` to your
account's rate limit. The marker/ledger/validation machinery is unchanged.

## Notes
- Outputs are **silver** — schema-valid ≠ signed off. Promotion to `reference_records/` stays a
  deliberate human step (cell-by-cell audit), exactly as for Rojas/Keskin.
- Cost: ~$3–8/paper on Opus-1M (the configured model). The ledger sums it; watch it.
- Don't bake API keys/creds into the image or commit them — they're bind-mounted/env-injected.
