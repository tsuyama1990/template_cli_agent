#!/bin/bash
set -e

# Wrapper for the Python contract generator
# Usage: ./scripts/ai_gen_contract.sh "User model with name and age"

if [ -z "$1" ]; then
    echo "Usage: $0 \"<description>\""
    exit 1
fi

uv run python scripts/ai_controller.py gen-contract "$1"
