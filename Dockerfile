# vaxtract — BYOK neoantigen cancer-vaccine extraction agent.
#
# Full image (Claude CLI + OpenAI-compat runtimes):
#   docker build --target full -t ghcr.io/mskgreenbaumlab/vaxtract:latest .
# Slim image (Grok/GPT/vLLM only — no Node/Claude CLI):
#   docker build --target compat -t ghcr.io/mskgreenbaumlab/vaxtract:compat .
#
# Run:
#   docker run --rm -e ANTHROPIC_API_KEY -v "$PWD/paper:/work/paper" \
#       ghcr.io/mskgreenbaumlab/vaxtract:latest /work/paper /work/paper/out.json
#   docker run --rm -e XAI_API_KEY -v "$PWD/paper:/work/paper" \
#       ghcr.io/mskgreenbaumlab/vaxtract:compat \
#       --profile grok /work/paper /work/paper/out.json

FROM python:3.11-slim AS base
WORKDIR /app
COPY pyproject.toml README.md LICENSE /app/
COPY vaxtract /app/vaxtract

# Slim: OpenAI-compatible profiles only (Grok, GPT, Inkling, NIM/vLLM).
FROM base AS compat
RUN pip install --no-cache-dir ".[compat_figures]"
WORKDIR /work
ENTRYPOINT ["vaxtract"]

# Default: Claude Agent SDK + figures + OpenAI-compat. Bundles Node + `claude`.
FROM base AS full
RUN apt-get update \
 && apt-get install -y --no-install-recommends curl ca-certificates gnupg \
 && curl -fsSL https://deb.nodesource.com/setup_20.x | bash - \
 && apt-get install -y --no-install-recommends nodejs \
 && npm install -g @anthropic-ai/claude-code \
 && apt-get purge -y curl gnupg && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir ".[agent,figures,compat]"
WORKDIR /work
ENTRYPOINT ["vaxtract"]
