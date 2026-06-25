# vaxtract — BYOK neoantigen cancer-vaccine extraction agent.
#
# Build:  docker build -t vaxtract .
# Run:    docker run --rm -e ANTHROPIC_API_KEY \
#               -v "$PWD/paper:/work/paper" \
#               vaxtract /work/paper /work/paper/out.json
#
# The image bundles Node + the Claude Code CLI because the Claude Agent SDK shells
# out to the `claude` binary. Bring your own ANTHROPIC_API_KEY (nothing is hosted).
FROM python:3.11-slim

# Node 20 + the Claude Code CLI (the SDK's subprocess transport).
RUN apt-get update \
 && apt-get install -y --no-install-recommends curl ca-certificates gnupg \
 && curl -fsSL https://deb.nodesource.com/setup_20.x | bash - \
 && apt-get install -y --no-install-recommends nodejs \
 && npm install -g @anthropic-ai/claude-code \
 && apt-get purge -y curl gnupg && apt-get autoremove -y \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app
# Copy only what the build needs (see .dockerignore); keeps the image lean.
COPY pyproject.toml README_PACKAGE.md LICENSE /app/
COPY vaxtract /app/vaxtract
# This image runs the agent, so it needs the [agent] extra (SDK + readers);
# [figures] adds PDF figure/image reading. Core install alone is schema-only.
RUN pip install --no-cache-dir ".[agent,figures]"

# Papers are mounted here; pass <paper_dir> [out.json] as args.
WORKDIR /work
ENTRYPOINT ["vaxtract"]
