# Model-agnostic extraction runtime

The extraction **record** (schema, `agent_core`, prompt) was already SDK-free. The
**loop** now is too: two runtimes, many named profiles.

| Runtime | How you reach it | Models |
|---|---|---|
| `claude_sdk` | Claude Agent SDK (default) | Opus / Sonnet / Fable via `--model` |
| `openai_compat` | Chat Completions + tools, `base_url` + API key | Grok, GPT, Inkling, Nemotron, Qwen, any OpenAI-compatible server |

Non-Claude profiles are **experimental**. Schema-valid JSON is not enough to
call a profile supported.

## CLI

```bash
pip install 'vaxtract[agent]'     # Claude
pip install 'vaxtract[compat]'    # openai Python client (Grok / GPT / vLLM / Together)

vaxtract ./paper out.json                              # Claude default
vaxtract --subscription --model claude-opus-5[1m] ./paper out.json
vaxtract --profile grok ./paper out.json               # XAI_API_KEY
vaxtract --profile gpt --model gpt-5.4 ./paper out.json
vaxtract --profile inkling ./paper out.json            # TOGETHER_API_KEY
vaxtract --backend openai_compat --base-url http://127.0.0.1:8000/v1 \
         --model nvidia/NVIDIA-Nemotron-3-Ultra-550B-A55B ./paper out.json
```

`--subscription` is Claude-only. Combining it with `--profile grok` (or any
`openai_compat` runtime) is a hard error.

Profiles live in `vaxtract/profiles.yaml`. `vision: false` makes `read_figure`
return an error instead of inventing bar heights.

## Eval

Score an extraction against the audited `reference_records/` (not a live model).
Live `--profile grok` runs in 2026-08-23 matched Claude on peptide/epitope identity
for Cafri and Keskin but missed the evidence-recall floor. `grok` stays experimental.

## Docker

```bash
# full: Claude CLI + OpenAI-compat
docker pull ghcr.io/mskgreenbaumlab/vaxtract:latest
# slim: Grok/GPT only
docker pull ghcr.io/mskgreenbaumlab/vaxtract:compat
```

Images are published from GitHub Releases (`docker.yml`), the same event that
publishes the PyPI wheel.
