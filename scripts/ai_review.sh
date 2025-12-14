#!/bin/bash
set -e

# Wrapper for the Python strict review controller
uv run python scripts/ai_controller.py strict-review
