"""This module contains utility functions for the AC-CDD application."""

from pathlib import Path


def detect_package_dir() -> str:
    """
    Detects the main package directory under dev_src/.
    Looks for the first directory containing __init__.py.
    """
    src_path = Path("dev_src")
    if src_path.exists():
        for p in src_path.iterdir():
            if p.is_dir() and (p / "__init__.py").exists():
                return str(p)
    return "dev_src/ac_cdd_core"


PROMPT_FILENAME_MAP = {
    "auditor.md": "AUDITOR_INSTRUCTION.md",
    "coder.md": "CODER_INSTRUCTION.md",
    "architect.md": "ARCHITECT_INSTRUCTION.md",
}


def read_prompt(filename: str, default: str) -> str:
    """
    Reads a prompt file from various locations, with a fallback.
    """
    target_filename = PROMPT_FILENAME_MAP.get(filename, filename)

    user_prompt_mapped = Path("dev_documents/system_prompts") / target_filename
    if user_prompt_mapped.exists():
        return user_prompt_mapped.read_text(encoding="utf-8").strip()

    user_prompt_direct = Path("dev_documents/system_prompts") / filename
    if user_prompt_direct.exists():
        return user_prompt_direct.read_text(encoding="utf-8").strip()

    system_prompt = Path("dev_src/ac_cdd_core/prompts") / filename
    if system_prompt.exists():
        return system_prompt.read_text(encoding="utf-8").strip()

    return default
