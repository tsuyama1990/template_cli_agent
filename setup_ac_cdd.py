import os
import shutil
import sys
from pathlib import Path
import subprocess

# 定数定義
BASE_DIR = Path(__file__).resolve().parent

def create_structure():
    """
    templates ディレクトリからプロジェクト構造を作成する
    """
    print("Setting up AC-CDD Environment from templates...")

    # Template source
    templates_root = BASE_DIR / "templates"
    
    if not templates_root.exists():
        print(f"Error: templates directory not found at {templates_root}")
        sys.exit(1)

    # 1. Directories to create (ensure they exist)
    dirs = [
        "src/contracts",
        "src/core",
        "tests/property",
        "tests/e2e",
        "documents/templates",
        "scripts",
        ".github/workflows",
    ]
    for d in dirs:
        Path(d).mkdir(parents=True, exist_ok=True)
        print(f"Created directory: {d}")

    # 2. Files Mapping: Source (relative to templates/) -> Dest (relative to root)
    # Note: simple copy logic
    
    # Files to copy
    # (Source in templates, Destination in project)
    files_to_copy = [
        ("pyproject.toml", "pyproject.toml"),
        ("manage.py", "manage.py"),
        ("env.example", ".env.example"),
        ("scripts/ai_orchestrator.py", "scripts/ai_orchestrator.py"),
        ("scripts/utils.py", "scripts/utils.py"),
        ("src/main.py", "src/main.py"),
        ("contracts/schema.py", "documents/templates/schema_template.py"),
        ("documents/SPEC_TEMPLATE.md", "documents/templates/SPEC_TEMPLATE.md"),
        ("documents/UAT_TEMPLATE.md", "documents/templates/UAT_TEMPLATE.md"),
        ("github/workflows/ci.yml", ".github/workflows/ci.yml"),
        # We also copy ac_cdd.toml if it exists in templates, to allow easy config
    ]
    
    if (templates_root / "ac_cdd.toml").exists():
        files_to_copy.append(("ac_cdd.toml", "ac_cdd.toml"))

    for src_rel, dst_rel in files_to_copy:
        src = templates_root / src_rel
        dst = Path(dst_rel)
        if src.exists():
            shutil.copy(src, dst)
            print(f"Created file: {dst_rel}")
        else:
            print(f"Warning: Template file missing: {src_rel}")

    # 3. Init files
    Path("src/contracts/__init__.py").touch()
    Path("scripts/__init__.py").touch()

    # 4. Make manage.py executable
    try:
        mode = os.stat("manage.py").st_mode
        os.chmod("manage.py", mode | 0o755)
    except Exception as e:
        print(f"Warning: Could not make manage.py executable: {e}")
        
    print("Files copied successfully.")

if __name__ == "__main__":
    create_structure()
    print("Setup Complete! Run 'make install' (or 'uv sync') and then './manage.py init' to start.")
