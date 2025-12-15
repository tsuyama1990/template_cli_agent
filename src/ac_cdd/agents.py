import os
from dataclasses import dataclass
from pathlib import Path

from pydantic_ai import Agent, RunContext
from pydantic_ai.models.gemini import GeminiModel

from .domain_models import CyclePlan, AuditResult, UatAnalysis, FileArtifact

# Configure the model
# Ensure GEMINI_API_KEY is available in environment for pydantic-ai
# Using gemini-2.0-flash as per latest available models
model = GeminiModel('gemini-2.0-flash')

@dataclass
class AgentDeps:
    documents_dir: Path
    cycle_id: str | None = None

# Personas
ARCHITECT_PERSONA = (
    "あなたはシニア・ソフトウェアアーキテクトです。"
    "要件定義書に基づき、堅牢でスケーラブルな設計仕様を策定してください。"
)
CODER_PERSONA = (
    "あなたはJules、熟練したPythonエンジニアです。"
    "与えられた仕様と契約(Contract)に基づき、高品質なコードを実装してください。"
    "必ず `<thought>` タグで思考プロセスを出力してからコードを書いてください。"
)
TESTER_PERSONA = (
    "あなたはQAエンジニアです。"
    "実装の詳細には立ち入らず、契約(Contract)のみに基づいて、エッジケースを網羅するテストケースを作成してください。"
)
AUDITOR_PERSONA = (
    "あなたは世界一厳格なコード監査人(Gemini)です。"
    "Pydantic契約違反、セキュリティ、設計原則の観点からコードを徹底的にレビューしてください。"
)
QA_ANALYST_PERSONA = (
    "あなたはQAマネージャーです。"
    "テストログを分析し、要件に対する適合度と挙動の考察をMarkdownで報告してください。"
)

# Helper for dynamic context
def _get_common_context(deps: AgentDeps) -> str:
    context = []

    # Constitution
    all_spec = deps.documents_dir / "ALL_SPEC.md"
    if all_spec.exists():
        context.append(f"=== Project Constitution (ALL_SPEC.md) ===\n{all_spec.read_text(encoding='utf-8')}")

    # Conventions
    conventions = deps.documents_dir / "conventions.md"
    if conventions.exists():
        context.append(f"=== Conventions ===\n{conventions.read_text(encoding='utf-8')}")

    return "\n\n".join(context)

# --- Agents ---

# Planner
planner_agent = Agent(
    model,
    output_type=CyclePlan,
    system_prompt=ARCHITECT_PERSONA,
    deps_type=AgentDeps,
)

@planner_agent.system_prompt
def planner_context(ctx: RunContext[AgentDeps]) -> str:
    return _get_common_context(ctx.deps)

# Coder
coder_agent = Agent(
    model,
    output_type=list[FileArtifact],
    system_prompt=CODER_PERSONA,
    deps_type=AgentDeps,
)

@coder_agent.system_prompt
def coder_context(ctx: RunContext[AgentDeps]) -> str:
    return _get_common_context(ctx.deps)

# Tester (returns list of files to support multiple test files if needed)
tester_agent = Agent(
    model,
    output_type=list[FileArtifact],
    system_prompt=TESTER_PERSONA,
    deps_type=AgentDeps,
)

@tester_agent.system_prompt
def tester_context(ctx: RunContext[AgentDeps]) -> str:
    return _get_common_context(ctx.deps)

# Auditor
auditor_agent = Agent(
    model,
    output_type=AuditResult,
    system_prompt=AUDITOR_PERSONA,
    deps_type=AgentDeps,
)

@auditor_agent.system_prompt
def auditor_context(ctx: RunContext[AgentDeps]) -> str:
    return _get_common_context(ctx.deps)

# QA Analyst
qa_agent = Agent(
    model,
    output_type=UatAnalysis,
    system_prompt=QA_ANALYST_PERSONA,
    deps_type=AgentDeps,
)

@qa_agent.system_prompt
def qa_context(ctx: RunContext[AgentDeps]) -> str:
    return _get_common_context(ctx.deps)
