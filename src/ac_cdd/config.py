import tomllib
from pathlib import Path

from pydantic import ConfigDict
from pydantic_settings import BaseSettings, SettingsConfigDict


class PathsConfig(BaseSettings):
    model_config = ConfigDict(extra="ignore")
    documents_dir: str = "dev_documents"
    contracts_dir: str = "src/ac_cdd/contracts"
    sessions_dir: str = ".jules/sessions"
    src: str = "src"
    tests: str = "tests"
    templates: str = "dev_documents/templates"

class ToolsConfig(BaseSettings):
    model_config = ConfigDict(extra="ignore")
    # Deprecated jules_cmd/gemini_cmd/jules_base_url but keeping struct if needed for compat
    # or we can remove them. I'll remove deprecated ones to be clean.
    gh_cmd: str = "gh"
    audit_cmd: str = "bandit"
    uv_cmd: str = "uv"
    mypy_cmd: str = "mypy"

class AgentsConfig(BaseSettings):
    # Agents are now defined in code (agents.py) but we might want to keep config for something?
    # The requirement said "config.py から読み込むのではなく...".
    # So we can remove the prompts from here.
    model_config = ConfigDict(extra="ignore")
    # Keeping empty or removing entirely?
    # The Settings class instantiates it, so let's keep it empty/flexible or remove field.
    pass

class PromptsConfig(BaseSettings):
    model_config = ConfigDict(extra="ignore")
    property_test_template: str = "実装は見ず、このPydanticスキーマ (contracts/) の制約が正しく機能するかを検証するHypothesisテストを作成せよ。出力先は tests/property/test_cycle{cycle_id}.py"

class Settings(BaseSettings):
    MAX_RETRIES: int = 10

    paths: PathsConfig = PathsConfig()
    tools: ToolsConfig = ToolsConfig()
    # agents: AgentsConfig = AgentsConfig() # Removing as unused
    prompts: PromptsConfig = PromptsConfig()

    model_config = SettingsConfigDict(
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore"
    )

    @classmethod
    def load_from_toml(cls, toml_path: str = "ac_cdd.toml") -> "Settings":
        settings = cls()
        path = Path(toml_path)
        if path.exists():
            with open(path, "rb") as f:
                data = tomllib.load(f)
            if "paths" in data:
                settings.paths = PathsConfig(**data["paths"])
            if "tools" in data:
                settings.tools = ToolsConfig(**data["tools"])
            # agents removed
            if "prompts" in data:
                settings.prompts = PromptsConfig(**data["prompts"])
        return settings

settings = Settings.load_from_toml()
