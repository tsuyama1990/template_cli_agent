import json
import os
import shutil
import subprocess
import time
import asyncio
from pathlib import Path

from .config import settings
from .tools import ToolNotFoundError, ToolWrapper
from .utils import logger
from .agents import (
    AgentDeps,
    planner_agent,
    coder_agent,
    tester_agent,
    auditor_agent,
    qa_agent
)
from .domain_models import CyclePlan, AuditResult, UatAnalysis, FileArtifact


class CycleOrchestrator:
    """
    AC-CDD サイクルの自動化を管理するクラス。
    実装、テスト、監査、UATの各フェーズをオーケストレーションする。
    """

    def __init__(self, cycle_id: str, dry_run: bool = False, auto_next: bool = False) -> None:
        self.cycle_id = cycle_id
        self.dry_run = dry_run
        self.auto_next = auto_next

        # Use paths from config
        self.documents_dir = Path(settings.paths.documents_dir)
        self.contracts_dir = Path(settings.paths.contracts_dir)

        self.cycle_dir = self.documents_dir / f"CYCLE{cycle_id}"
        self.audit_log_path = self.cycle_dir / "AUDIT_LOG.md"

        if not self.cycle_dir.exists():
            # In new cycle, it might not exist yet, that's fine if we are creating it
            pass

        # Initialize tools
        try:
            self.gh = ToolWrapper(settings.tools.gh_cmd)
            self.audit_tool = ToolWrapper(settings.tools.audit_cmd) # bandit
            self.uv = ToolWrapper(settings.tools.uv_cmd)
            self.mypy = ToolWrapper(settings.tools.mypy_cmd)
        except ToolNotFoundError as e:
            if not self.dry_run:
                raise
            logger.warning(f"[DRY-RUN] Tool missing: {e}. Proceeding anyway.")

    async def execute_all(self, progress_task=None, progress_obj=None) -> None:
        """全フェーズを実行"""
        # Note: steps should be async aware
        steps = [
            ("Planning Phase", self.plan_cycle),
            ("Aligning Contracts", self.align_contracts),
            ("Generating Property Tests", self.generate_property_tests),
            ("Implementation Loop", self.run_implementation_loop),
            ("UAT Phase", self.run_uat_phase),
            ("Finalizing Cycle", self.finalize_cycle),
        ]

        for name, func in steps:
            if progress_obj:
                progress_obj.update(progress_task, description=f"[cyan]{name}...")
            logger.info(f"Starting Phase: {name}")
            if asyncio.iscoroutinefunction(func):
                await func()
            else:
                func()
            logger.info(f"Completed Phase: {name}")

    async def plan_cycle(self) -> None:
        """
        Phase 0: 計画策定 (Planning Phase)
        JulesがALL_SPEC.mdとCYCLE_PLANNING_PROMPT.mdを読み込み、
        自動的に実装計画 (SPEC.md, schema.py, UAT.md) を策定・配置する。
        """
        if self.dry_run:
            logger.info("[DRY-RUN] Planning Cycle... (Mocking plan generation)")
            return

        # 1. Prepare inputs
        planning_prompt_path = Path(settings.paths.templates) / "CYCLE_PLANNING_PROMPT.md"
        base_prompt = planning_prompt_path.read_text() if planning_prompt_path.exists() else ""

        user_task = base_prompt + f"\n\nFocus specifically on generating artifacts for CYCLE{self.cycle_id}."

        deps = AgentDeps(documents_dir=self.documents_dir, cycle_id=self.cycle_id)

        # 2. Call Agent
        logger.info(f"Generating Plan for CYCLE{self.cycle_id}...")

        result = await planner_agent.run(user_task, deps=deps)
        plan: CyclePlan = result.data

        # 3. Save Artifacts
        self.cycle_dir.mkdir(parents=True, exist_ok=True)

        # Helper to save artifact
        def save_artifact(artifact: FileArtifact):
            # Enforce saving in cycle dir
            fname = Path(artifact.path).name
            target = self.cycle_dir / fname
            target.write_text(artifact.content)
            logger.info(f"Saved {target}")

        save_artifact(plan.spec_file)
        save_artifact(plan.schema_file)
        save_artifact(plan.uat_file)

        # Save thought process
        (self.cycle_dir / "PLAN_THOUGHTS.md").write_text(plan.thought_process)

    def align_contracts(self) -> None:
        """
        Phase 3.1: 契約の整合確認とマージ
        src/ac_cdd/contracts/ に Cycleのschema.pyをマージする。
        """
        source_schema = self.cycle_dir / "schema.py"
        target_schema = self.contracts_dir / f"schema_cycle{self.cycle_id}.py"

        if not source_schema.exists():
            # Warning only if we are not strict or if previous step failed silently
            if not self.dry_run:
                 raise FileNotFoundError(f"{source_schema} not found.")
            else:
                 logger.warning(f"[DRY-RUN] Source schema {source_schema} missing.")
                 return

        if self.dry_run:
            logger.info(f"[DRY-RUN] Would copy {source_schema} to {target_schema}")
            return

        self.contracts_dir.mkdir(parents=True, exist_ok=True)

        # 既存ファイルがある場合、バックアップ
        if target_schema.exists():
            backup = target_schema.with_suffix(".py.bak")
            shutil.copy(target_schema, backup)
            logger.info(f"Backed up existing schema to {backup}")

        shutil.copy(source_schema, target_schema)

        # __init__.py の更新
        init_file = self.contracts_dir / "__init__.py"
        import_line = f"from .schema_cycle{self.cycle_id} import *"

        if init_file.exists():
            content = init_file.read_text()
            if import_line not in content:
                with open(init_file, "a") as f:
                    f.write(f"\n{import_line}\n")
        else:
            with open(init_file, "w") as f:
                f.write(f"{import_line}\n")

    async def generate_property_tests(self) -> None:
        """
        Phase 3.2: プロパティベーステストの生成
        """
        user_task = settings.prompts.property_test_template.format(cycle_id=self.cycle_id)
        deps = AgentDeps(documents_dir=self.documents_dir, cycle_id=self.cycle_id)

        if self.dry_run:
            logger.info(f"[DRY-RUN] Generating property tests: {user_task}")
            return

        # Call Agent (using tester agent)
        result = await tester_agent.run(user_task, deps=deps)
        artifacts: list[FileArtifact] = result.data

        for artifact in artifacts:
            # Save to tests/property/
            target_path = Path("tests/property") / Path(artifact.path).name
            target_path.parent.mkdir(parents=True, exist_ok=True)
            target_path.write_text(artifact.content)
            logger.info(f"Saved test artifact: {target_path}")

    async def run_implementation_loop(self) -> None:
        """
        Phase 3.3 & 3.4: 実装・CI・監査ループ
        """
        logger.info("Starting Implementation Phase")

        max_plan_retries = 3
        plan_attempt = 0

        while plan_attempt < max_plan_retries:
            plan_attempt += 1
            if plan_attempt > 1:
                logger.warning(f"Self-Healing: Re-planning cycle ({plan_attempt}/{max_plan_retries})...")
                await self._replan_cycle("Implementation loop failed repeatedly. Review SPEC/Schema/UAT.")

            # 1. Initial Implementation
            await self._trigger_implementation()

            # 2. Refinement Loop
            max_retries = settings.MAX_RETRIES
            attempt = 0

            logger.info("Entering Refinement Loop (Stable Audit Loop)")
            loop_success = False

            while attempt < max_retries:
                attempt += 1
                logger.info(f"Refinement Loop: Iteration {attempt}/{max_retries}")

                # 2.1 Test
                logger.info("Running Tests...")
                passed, logs = self._run_tests()

                if not passed:
                    logger.warning("Tests Failed. Triggering fix...")
                    fix_prompt = (
                        "Test Failed.\n"
                        f"Logs (truncated):\n{logs[-2000:]}\n"
                        "Please analyze the stack trace and fix the implementation in src/."
                    )
                    await self._trigger_fix(fix_prompt)
                    continue

                # 2.2 Audit
                logger.info("Tests Passed. Proceeding to Strict Audit...")
                audit_passed = await self.run_strict_audit()

                if audit_passed:
                    logger.info("Audit Passed (Clean)!")
                    loop_success = True
                    break
                else:
                    logger.warning("Audit Failed. Triggering fix...")
                    # Fix based on audit log
                    await self._trigger_fix(f"Audit failed. See {self.audit_log_path} for details.")
                    continue

            if loop_success:
                return

            logger.error("Implementation Loop Failed.")
            # Continue to outer loop (Self-Healing)

        raise Exception(f"Max retries reached in Self-Healing Plan Loop ({max_plan_retries} attempts).")

    async def _replan_cycle(self, feedback: str) -> None:
        """Re-runs planning with feedback."""
        logger.info("Triggering Re-Planning with feedback...")

        planning_prompt_path = Path(settings.paths.templates) / "CYCLE_PLANNING_PROMPT.md"
        base_prompt = planning_prompt_path.read_text() if planning_prompt_path.exists() else ""

        user_task = base_prompt + (
            f"\n\nCRITICAL UPDATE: The previous plan failed during implementation.\n"
            f"Feedback: {feedback}\n"
            f"Please revise the SPEC, Schema, and UAT for CYCLE{self.cycle_id}."
        )

        deps = AgentDeps(documents_dir=self.documents_dir, cycle_id=self.cycle_id)
        result = await planner_agent.run(user_task, deps=deps)
        plan: CyclePlan = result.data

        # Overwrite artifacts
        (self.cycle_dir / "SPEC.md").write_text(plan.spec_file.content)
        (self.cycle_dir / "schema.py").write_text(plan.schema_file.content)
        (self.cycle_dir / "UAT.md").write_text(plan.uat_file.content)

    async def _trigger_implementation(self) -> None:
        spec_path = f"{settings.paths.documents_dir}/CYCLE{self.cycle_id}/SPEC.md"
        description = (
            f"Implement requirements in {spec_path} "
            f"following schema in {settings.paths.contracts_dir}/"
        )

        deps = AgentDeps(documents_dir=self.documents_dir, cycle_id=self.cycle_id)

        if self.dry_run:
            logger.info(f"[DRY-RUN] Implementing feature: {description}")
            return

        result = await coder_agent.run(description, deps=deps)
        files: list[FileArtifact] = result.data

        for f in files:
            # Ensure path is relative to src or root?
            # The agent might return "src/foo/bar.py".
            # We trust the agent to return valid relative paths.
            p = Path(f.path)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(f.content)
            logger.info(f"Implemented: {p}")

    def _run_tests(self) -> tuple[bool, str]:
        """Runs tests locally using uv run pytest."""
        if self.dry_run:
            logger.info("[DRY-RUN] Running tests locally... (Mocking success)")
            return True, ""

        uv_path = shutil.which("uv")
        if not uv_path:
            raise ToolNotFoundError("uv not found")

        cmd = [uv_path, "run", "pytest"]
        try:
            result = subprocess.run( # noqa: S603
                cmd, capture_output=True, text=True, check=False
            )
            logs = result.stdout + "\n" + result.stderr
            return result.returncode == 0, logs
        except Exception as e:
            logger.error(f"Failed to run tests: {e}")
            return False, str(e)

    async def _trigger_fix(self, instructions: str) -> None:
        if self.dry_run:
            logger.info(f"[DRY-RUN] Fixing code: {instructions}")
            return

        deps = AgentDeps(documents_dir=self.documents_dir, cycle_id=self.cycle_id)
        # Use coder agent for fixes as well
        # We might want to pass message history?
        # For now, stateless fix request.
        result = await coder_agent.run(instructions, deps=deps)
        files: list[FileArtifact] = result.data

        for f in files:
            p = Path(f.path)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(f.content)
            logger.info(f"Fixed: {p}")

    async def run_strict_audit(self) -> bool:
        """
        Phase 3.4: 世界一厳格な監査
        """
        if self.dry_run:
            logger.info("[DRY-RUN] Audit... (Mocking approval)")
            return True

        logger.info("Running Static Analysis (Ruff, Mypy, Bandit)...")

        # 1. Ruff
        ruff_path = shutil.which("ruff")
        if ruff_path:
            try:
                subprocess.run([ruff_path, "check", "--fix", "src/", "tests/"], check=True, capture_output=True) # noqa: S603
                subprocess.run([ruff_path, "format", "src/", "tests/"], check=True, capture_output=True) # noqa: S603
            except subprocess.CalledProcessError as e:
                msg = f"Linting (Ruff) Failed:\n{e.stderr}"
                self._log_audit_failure([msg])
                return False

        # 2. Mypy
        try:
            self.mypy.run(["src/"])
        except Exception as e:
            self._log_audit_failure([f"Mypy Failed: {e}"])
            return False

        # 3. Bandit
        try:
            self.audit_tool.run(["-r", "src/", "-ll"])
        except Exception as e:
            self._log_audit_failure([f"Bandit Failed: {e}"])
            return False

        # 4. LLM Audit
        logger.info("Static checks passed. Proceeding to LLM Audit...")
        files_to_audit = self._get_filtered_files("src/")
        files_content = ""
        for fpath in files_to_audit:
            try:
                files_content += f"\n\n=== File: {fpath} ===\n{Path(fpath).read_text()}"
            except Exception:
                pass

        user_task = f"Audit the following code files:\n{files_content}"
        deps = AgentDeps(documents_dir=self.documents_dir, cycle_id=self.cycle_id)

        result = await auditor_agent.run(user_task, deps=deps)
        audit_res: AuditResult = result.data

        if audit_res.is_approved:
            return True
        else:
            self._log_audit_failure(audit_res.critical_issues + audit_res.suggestions)
            return False

    def _get_filtered_files(self, directory: str) -> list[str]:
        ignored_patterns = {"__pycache__", ".git", ".env", ".DS_Store", "*.pyc"}
        auditignore_path = Path(".auditignore")
        if auditignore_path.exists():
            for line in auditignore_path.read_text().splitlines():
                if line.strip() and not line.startswith("#"):
                    ignored_patterns.add(line.strip())

        import fnmatch
        files = []
        path = Path(directory)
        for p in path.rglob("*"):
            if p.is_file():
                is_ignored = False
                for pattern in ignored_patterns:
                    if fnmatch.fnmatch(p.name, pattern) or fnmatch.fnmatch(str(p), pattern) or pattern in str(p):
                        is_ignored = True
                        break
                if not is_ignored:
                    files.append(str(p))
        return files

    def _log_audit_failure(self, comments: list[str]) -> None:
        with open(self.audit_log_path, "a") as f:
            f.write(f"\n## Audit Failed at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            for c in comments:
                f.write(f"- {c}\n")

    async def run_uat_phase(self) -> None:
        """
        Phase 3.5: UATの生成と実行
        """
        if self.dry_run:
            logger.info("[DRY-RUN] UAT... (Mocking success)")
            return

        # 1. Generate UAT Code
        uat_path = f"{settings.paths.documents_dir}/CYCLE{self.cycle_id}/UAT.md"
        description = (
            f"Create Playwright tests in tests/e2e/ based on {uat_path}.\n"
            "REQUIREMENTS:\n"
            "- Use `unittest.mock`, `pytest-mock`, or `vcrpy` to mock ALL external connections.\n"
            "- Output valid Python code."
        )
        deps = AgentDeps(documents_dir=self.documents_dir, cycle_id=self.cycle_id)

        result = await tester_agent.run(description, deps=deps)
        files: list[FileArtifact] = result.data
        for f in files:
            p = Path(f.path)
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(f.content)

        # 2. Run Tests
        uv_path = shutil.which("uv")
        if not uv_path:
            raise ToolNotFoundError("uv not found")

        cmd = [uv_path, "run", "pytest", "tests/e2e/"]
        success = False
        logs = ""
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, check=False) # noqa: S603
            logs = res.stdout + "\n" + res.stderr
            success = (res.returncode == 0)
        except Exception as e:
            logs = str(e)

        # 3. Analyze
        await self._analyze_uat_results(logs, success)
        if not success:
            await self._trigger_fix(f"UAT Tests Failed. Logs:\n{logs[-2000:]}")
            raise Exception("UAT Phase Failed")

    async def _analyze_uat_results(self, logs: str, success: bool) -> None:
        logger.info("Analyzing UAT Results...")
        verdict = "PASS" if success else "FAIL"

        user_task = (
            f"Analyze the following pytest logs for UAT. Verdict: {verdict}\n"
            f"Logs (truncated): {logs[-10000:]}"
        )
        deps = AgentDeps(documents_dir=self.documents_dir, cycle_id=self.cycle_id)

        result = await qa_agent.run(user_task, deps=deps)
        analysis: UatAnalysis = result.data

        report_path = self.cycle_dir / "UAT_RESULT.md"
        with open(report_path, "w") as f:
            f.write(f"# UAT Result: {analysis.verdict}\n\n")
            f.write(f"## Summary\n{analysis.summary}\n\n")
            f.write(f"## Behavior Analysis\n{analysis.behavior_analysis}\n")
        logger.info(f"UAT Report saved to {report_path}")

    async def finalize_cycle(self) -> None:
        """
        Phase 4: 自動マージ
        """
        if self.dry_run:
            logger.info("[DRY-RUN] Merging PR via gh CLI...")
        else:
            args = ["pr", "merge", "--squash", "--delete-branch", "--admin"]
            self.gh.run(args)

        if self.auto_next:
            await self.prepare_next_cycle()

    async def prepare_next_cycle(self) -> None:
        """
        Auto-Next: Scaffolds and starts planning for the next cycle.
        """
        logger.info(f"Auto-Next enabled: Preparing next cycle after CYCLE{self.cycle_id}...")

        try:
            current_int = int(self.cycle_id)
            next_id = f"{current_int + 1:02d}"
        except ValueError:
            return

        next_cycle_dir = self.documents_dir / f"CYCLE{next_id}"
        if not next_cycle_dir.exists():
            next_cycle_dir.mkdir(parents=True)
            templates_dir = Path(settings.paths.templates) / "cycle"
            for item in ["SPEC.md", "UAT.md", "schema.py"]:
                src = templates_dir / item
                dst = next_cycle_dir / item
                if src.exists():
                    shutil.copy(src, dst)

        logger.info(f"Triggering Planning Phase for CYCLE{next_id}...")

        next_orchestrator = CycleOrchestrator(
            next_id, dry_run=self.dry_run, auto_next=self.auto_next
        )
        await next_orchestrator.plan_cycle()
        logger.info(f"CYCLE{next_id} Planning Complete.")
