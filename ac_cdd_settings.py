from typing import Literal

from dotenv import load_dotenv

load_dotenv()

from pydantic import BaseModel, Field, model_validator
from pydantic_settings import BaseSettings, SettingsConfigDict

from dev_src.ac_cdd_core.utils import detect_package_dir, read_prompt


class PathsConfig(BaseSettings):
    model_config = SettingsConfigDict(extra="ignore")
    documents_dir: str = "dev_documents"
    package_dir: str = Field(default_factory=detect_package_dir)
    contracts_dir: str = ""
    sessions_dir: str = ".jules/sessions"
    src: str = "src"
    tests: str = "tests"
    templates: str = "dev_documents/system_prompts"
    property_tests: str = "tests/property"
    unit_tests: str = "tests/unit"
    e2e_tests: str = "tests/e2e"
    prompts_dir: str = "dev_src/ac_cdd_core/prompts"

    @model_validator(mode="after")
    def _set_dependent_paths(self) -> "PathsConfig":
        if not self.contracts_dir:
            self.contracts_dir = f"{self.package_dir}/contracts"
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
        "pyproject.toml",
        "uv.lock",
        ".auditignore",
        "README.md",
        "ac_cdd_settings.py",
    ]
    install_cmd: str = "pip install --no-cache-dir ruff"
    test_cmd: list[str] = ["uv", "run", "pytest"]
    lint_check_cmd: list[str] = ["uv", "run", "ruff", "check", "--fix", "."]
    type_check_cmd: list[str] = ["uv", "run", "mypy", "src/"]
    security_check_cmd: list[str] = ["uv", "run", "bandit", "-r", "src/", "-ll"]


class AgentsConfig(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", extra="ignore")
    auditor_model: str = Field(
        default="openrouter/google/gemini-pro-1.5", validation_alias="SMART_MODEL"
    )
    qa_analyst_model: str = Field(
        default="openrouter/google/gemini-flash-1.5", validation_alias="FAST_MODEL"
    )


class ReviewerConfig(BaseSettings):
    model_config = SettingsConfigDict(extra="ignore")
    smart_model: str = Field(
        default="claude-3-5-sonnet",
        description="Model for editing code (Fixer)",
        validation_alias="SMART_MODEL",
    )
    fast_model: str = Field(
        default="openrouter/google/gemini-flash-1.5",
        description="Model for reading/auditing code",
        validation_alias="FAST_MODEL",
    )
    auditor: str = Field(
        default_factory=lambda: read_prompt("auditor.md", "DEFAULT_AUDITOR_PROMPT")
    )
    qa_analyst: str = Field(
        default_factory=lambda: read_prompt("UAT_DESIGN.md", "DEFAULT_QA_ANALYST_PROMPT")
    )


class ACCDDSettings(BaseSettings):
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

    class SessionConfig(BaseModel):
        session_id: str | None = None
        integration_branch_prefix: str = "dev"
        auto_merge_to_integration: bool = True
        final_merge_strategy: Literal["merge", "squash", "rebase"] = "squash"
        auto_delete_session_branches: bool = True

    session: SessionConfig = Field(default_factory=SessionConfig)
    paths: PathsConfig = PathsConfig()
    jules: JulesConfig = JulesConfig()
    tools: ToolsConfig = ToolsConfig()
    sandbox: SandboxConfig = SandboxConfig()
    agents: AgentsConfig = AgentsConfig()
    reviewer: ReviewerConfig = ReviewerConfig()

    @property
    def current_session_id(self) -> str:
        if self.session.session_id:
            return self.session.session_id
        from datetime import datetime
        now = datetime.now()
        return f"session-{now.strftime('%Y%m%d-%H%M%S-%f')[:20]}"

    @property
    def integration_branch(self) -> str:
        return f"{self.session.integration_branch_prefix}/{self.current_session_id}/integration"

    model_config = SettingsConfigDict(env_file=".env", env_file_encoding="utf-8", extra="allow")


config = ACCDDSettings()
