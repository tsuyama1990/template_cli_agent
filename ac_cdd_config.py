from pathlib import Path
import os

from dotenv import load_dotenv
from pydantic import Field, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

load_dotenv()  # Explicitly load .env into os.environ


def _detect_project_root() -> Path:
    """Detects the project root by looking for pyproject.toml."""
    current_path = Path.cwd()
    while current_path != current_path.parent:
        if (current_path / "pyproject.toml").exists():
            return current_path
        current_path = current_path.parent
    return Path.cwd()  # Fallback to current dir


def _detect_package_dir(project_root: Path) -> Path:
    """Detects the main package directory under dev_src/."""
    src_path = project_root / "dev_src"
    if src_path.exists():
        for p in src_path.iterdir():
            if p.is_dir() and (p / "__init__.py").exists():
                return p
    return src_path / "ac_cdd_core"


def _read_prompt(prompts_dir: Path, filename: str, default: str) -> str:
    p = prompts_dir / filename
    if p.exists():
        return p.read_text(encoding="utf-8").strip()
    return default


class PathsConfig(BaseSettings):
    model_config = SettingsConfigDict(extra="ignore")
    project_root: Path = Field(default_factory=_detect_project_root)
    package_dir: Path = Field(default=Path())  # Will be set in validator
    documents_dir: Path = Field(default=Path())
    contracts_dir: Path = Field(default=Path())
    sessions_dir: Path = Field(default=Path())
    src: Path = Field(default=Path())
    tests: Path = Field(default=Path())
    templates: Path = Field(default=Path())
    property_tests: Path = Field(default=Path())
    unit_tests: Path = Field(default=Path())
    e2e_tests: Path = Field(default=Path())
    prompts_dir: Path = Field(default=Path())

    @model_validator(mode="after")
    def _set_dependent_paths(self) -> "PathsConfig":
        if not self.package_dir.is_dir():
            self.package_dir = _detect_package_dir(self.project_root)
        if not self.documents_dir.is_dir():
            self.documents_dir = self.project_root / "dev_documents"
        if not self.contracts_dir.is_dir():
            self.contracts_dir = self.package_dir / "contracts"
        if not self.sessions_dir.is_dir():
            self.sessions_dir = self.project_root / ".jules/sessions"
        if not self.src.is_dir():
            self.src = self.project_root / "src"
        if not self.tests.is_dir():
            self.tests = self.project_root / "tests"
        if not self.templates.is_dir():
            self.templates = self.documents_dir / "templates"
        if not self.property_tests.is_dir():
            self.property_tests = self.tests / "property"
        if not self.unit_tests.is_dir():
            self.unit_tests = self.tests / "unit"
        if not self.e2e_tests.is_dir():
            self.e2e_tests = self.tests / "e2e"
        if not self.prompts_dir.is_dir():
            self.prompts_dir = self.package_dir / "prompts"
        return self


class JulesConfig(BaseSettings):
    model_config = SettingsConfigDict(extra="ignore")
    executable: str = "jules"
    timeout_seconds: int = 7200
    polling_interval_seconds: int = 120


class ToolsConfig(BaseSettings):
    model_config = SettingsConfigDict(extra="ignore")
    jules_cmd: str = "jules"
    gh_cmd: str = "gh"
    audit_cmd: str = "bandit"
    uv_cmd: str = "uv"
    mypy_cmd: str = "mypy"
    ruff_cmd: str = "ruff"
    gemini_cmd: str = "gemini"
    jules_base_url: str = "https://jules.googleapis.com/v1alpha"


class SandboxConfig(BaseSettings):
    model_config = SettingsConfigDict(extra="ignore")
    template: str | None = None
    timeout: int = 3600
    cwd: str = "/home/user/project"
    dirs_to_sync: list[str] = ["src", "tests", "contracts", "dev_documents", "dev_src"]
    files_to_sync: list[str] = [
        "pyproject.toml", "uv.lock", ".auditignore", "README.md", "ac_cdd_config.py"
    ]
    install_cmd: str = "pip install --no-cache-dir ruff aider-chat"
    test_cmd: list[str] = ["uv", "run", "pytest"]
    lint_check_cmd: list[str] = ["uv", "run", "ruff", "check", "--fix", "."]
    type_check_cmd: list[str] = ["uv", "run", "mypy", "src/"]
    security_check_cmd: list[str] = ["uv", "run", "bandit", "-r", "src/", "-ll"]


class AgentsConfig(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", extra="ignore")
    auditor_model: str = Field(default="openrouter/google/gemini-pro-1.5", validation_alias="SMART_MODEL")
    qa_analyst_model: str = Field(default="openrouter/google/gemini-flash-1.5", validation_alias="FAST_MODEL")


class AiderConfig(BaseSettings):
    model_config = SettingsConfigDict(extra="ignore")
    smart_model: str = Field(default="claude-3-5-sonnet", description="Model for editing code (Fixer)", validation_alias="SMART_MODEL")
    fast_model: str = Field(default="openrouter/google/gemini-flash-1.5", description="Model for reading/auditing code", validation_alias="FAST_MODEL")
    auditor: str = ""
    qa_analyst: str = ""


class PromptsConfig(BaseSettings):
    model_config = SettingsConfigDict(extra="ignore")
    property_test_template: str = ""
    structurer: str = ""


class Settings(BaseSettings):
    JULES_API_KEY: str | None = None
    OPENROUTER_API_KEY: str | None = None
    MAX_RETRIES: int = 10
    DUMMY_CYCLE_ID: str = "00"
    E2B_API_KEY: str | None = None
    GCP_PROJECT_ID: str | None = None
    GCP_REGION: str = "us-central1"
    NUM_AUDITORS: int = 3
    REVIEWS_PER_AUDITOR: int = 2
    MAX_ITERATIONS: int = 3

    paths: PathsConfig = Field(default_factory=PathsConfig)
    jules: JulesConfig = Field(default_factory=JulesConfig)
    tools: ToolsConfig = Field(default_factory=ToolsConfig)
    sandbox: SandboxConfig = Field(default_factory=SandboxConfig)
    agents: AgentsConfig = Field(default_factory=AgentsConfig)
    aider: AiderConfig = Field(default_factory=AiderConfig)
    prompts: PromptsConfig = Field(default_factory=PromptsConfig)

    model_config = SettingsConfigDict(env_file=".env", env_file_encoding="utf-8", extra="ignore")

    @model_validator(mode="after")
    def _load_prompts(self) -> "Settings":
        prompts_dir = self.paths.prompts_dir
        self.aider.auditor = _read_prompt(prompts_dir, "auditor.md", "DEFAULT_AUDITOR_PROMPT")
        self.aider.qa_analyst = _read_prompt(prompts_dir, "qa_analyst.md", "DEFAULT_QA_ANALYST_PROMPT")
        self.prompts.property_test_template = _read_prompt(prompts_dir, "property_test_template.md", "DEFAULT_TEST_PROMPT")
        self.prompts.structurer = str(prompts_dir / "structurer.md")
        return self


config = Settings()
