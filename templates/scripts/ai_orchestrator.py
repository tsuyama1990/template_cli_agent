import json
import shutil
import subprocess
import time
from pathlib import Path
try:
    import tomllib
except ImportError:
    # Fallback/Assumption: Python 3.11+
    try:
        import tomli as tomllib
    except ImportError:
        # Minimal TOML parser or fail if dep missing (should have been installed via project dep)
        # But this script runs inside the env, so dependencies should be there.
        # pyproject.toml adds 'tomli' for python < 3.11.
        import tomllib

from .utils import logger, run_command


class CycleOrchestrator:
    """
    AC-CDD サイクルの自動化を管理するクラス。
    実装、テスト、監査、UATの各フェーズをオーケストレーションする。
    """

    def __init__(self, cycle_id: str, dry_run: bool = False):
        self.cycle_id = cycle_id
        self.dry_run = dry_run
        self.cycle_dir = Path(f"documents/CYCLE{cycle_id}")
        self.contracts_dir = Path("src/contracts")
        self.audit_log_path = self.cycle_dir / "AUDIT_LOG.md"

        self.config = self._load_config()
        # Default keys
        self.ai_tool_command = self.config.get("ac_cdd", {}).get("ai_tool_command", ["jules", "remote", "new"])
        self.audit_tool_command = self.config.get("ac_cdd", {}).get("audit_tool_command", ["gemini", "--json"])

        if not self.cycle_dir.exists():
            raise ValueError(f"Cycle directory {self.cycle_dir} does not exist.")

    def _load_config(self):
        config_path = Path("ac_cdd.toml")
        if config_path.exists():
            try:
                with open(config_path, "rb") as f:
                    return tomllib.load(f)
            except Exception as e:
                logger.warning(f"Failed to load ac_cdd.toml: {e}. Using defaults.")
        return {}

    def execute_all(self, progress_task=None, progress_obj=None):
        """全フェーズを実行"""
        steps = [
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
            func()
            logger.info(f"Completed Phase: {name}")

    def align_contracts(self):
        """
        Phase 3.1: 契約の整合確認とマージ
        documents/CYCLE{id}/schema.py を src/contracts/ にマージする。
        """
        source_schema = self.cycle_dir / "schema.py"
        target_schema = self.contracts_dir / f"schema_cycle{self.cycle_id}.py"

        if not source_schema.exists():
            raise FileNotFoundError(f"{source_schema} not found.")

        # Pydanticモデルとして有効かチェック (Dry runでも行う)
        # 実際にはここで import してチェックしたいが、簡単のためファイルコピーとみなす

        if self.dry_run:
            logger.info(f"[DRY-RUN] Would copy {source_schema} to {target_schema}")
            return

        # ディレクトリ確認
        self.contracts_dir.mkdir(parents=True, exist_ok=True)

        # 既存ファイルがある場合、バックアップ
        if target_schema.exists():
            backup = target_schema.with_suffix(".py.bak")
            shutil.copy(target_schema, backup)
            logger.info(f"Backed up existing schema to {backup}")

        shutil.copy(source_schema, target_schema)

        # __init__.py の更新
        # 名前空間の衝突を避けるため、サイクルIDを含めたエイリアスでインポートするか、
        # あるいは「最新のサイクルが正義」として上書きを許容するか。
        # ここではシンプルに追記するが、実際の運用では古いサイクルのクラス名変更などを検討すべき。
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

    def generate_property_tests(self):
        """
        Phase 3.2: プロパティベーステストの生成
        Julesに契約のみを見せてテストを書かせる。
        """
        prompt = (
            "実装は見ず、このPydanticスキーマ (contracts/) の制約が正しく機能するかを検証する "
            "Hypothesisテストを作成せよ。出力先は tests/property/test_cycle{cycle_id}.py"
        )

        if self.dry_run:
            logger.info(f"[DRY-RUN] calling ai tool with prompt: {prompt}")
            return

        # contractsディレクトリの中身をコンテキストとして渡す
        cmd = self.ai_tool_command + [
            "--title", f"Generate Property Tests for Cycle {self.cycle_id}",
            "--description", prompt
            # コンテキスト指定が必要ならここに追加
        ]
        run_command(cmd)

    def run_implementation_loop(self):
        """
        Phase 3.3 & 3.4: 実装・CI・監査ループ
        """
        max_retries = 10
        attempt = 0

        while attempt < max_retries:
            attempt += 1
            logger.info(f"Implementation Loop: Attempt {attempt}/{max_retries}")

            # 1. Implement
            self._trigger_implementation()

            # 2. CI Check
            if not self._wait_for_ci():
                logger.warning("CI Failed. Triggering fix...")
                self._trigger_fix("CI Check failed. Please fix the implementation based on logs.")
                continue

            # 3. Strict Audit
            audit_result = self.run_strict_audit()
            if audit_result is True:
                logger.info("Audit Passed!")
                return
            else:
                logger.warning("Audit Failed. Triggering fix...")
                # Audit結果は AUDIT_LOG.md にある想定
                self._trigger_fix(f"Audit failed. See {self.audit_log_path} for details.")

        raise Exception("Max retries reached in Implementation Loop.")

    def _trigger_implementation(self):
        if self.dry_run:
            logger.info("[DRY-RUN] AI Implementing feature from SPEC.md and schema.py...")
            return

        cmd = self.ai_tool_command + [
            "--title",
            f"Implement Cycle {self.cycle_id}",
            "--description",
            f"Implement requirements in documents/CYCLE{self.cycle_id}/SPEC.md following "
            "schema in src/contracts/",
        ]
        run_command(cmd)

    def _wait_for_ci(self) -> bool:
        """GitHub Actionsの完了を待つ"""
        if self.dry_run:
            logger.info("[DRY-RUN] Waiting for CI... (Mocking success)")
            return True

        logger.info("Watching CI status...")
        # gh run watch を使う。失敗したら exit code が非0になる
        try:
            # 最新のRunを監視
            run_command(["gh", "run", "watch", "--exit-status"])
            return True
        except subprocess.CalledProcessError:
            return False

    def _trigger_fix(self, instructions: str):
        if self.dry_run:
            logger.info(f"[DRY-RUN] AI fixing code: {instructions}")
            return

        # Using configurable command. Appending instructions.
        # This implies a new session unless the tool supports state persistence or we are using a 'reply' tool.
        # For simplicity in this template, we assume 'new' works or user configures a tool that handles context.
        cmd = self.ai_tool_command + [
            "--title", f"Fix Cycle {self.cycle_id}",
            "--description", f"Fix request: {instructions}"
        ]
        run_command(cmd)

    def run_strict_audit(self) -> bool:
        """
        Phase 3.4: 世界一厳格な監査
        Gemini CLI (JSON mode) を呼び出す。
        """
        prompt = (
            "あなたは世界一厳格なコード監査人です。"
            "Pydantic契約違反、セキュリティ、設計原則の観点からコードをレビューしてください。"
            "合格なら `{\"approved\": true}`、"
            "不合格なら `{\"approved\": false, \"comments\": [...]}` をJSONで返してください。"
        )

        if self.dry_run:
            logger.info("[DRY-RUN] AI Auditing... (Mocking approval)")
            return True

        # srcディレクトリ全体を監査対象とする
        # gemini CLIの仕様に合わせて調整
        cmd = self.audit_tool_command + [
            prompt,
            "src/"    # 対象ディレクトリ
        ]

        try:
            # 実際には run_command で出力をキャプチャしてパースする
            # ここでは run_command はストリーミング出力用なので、subprocess.run を直接使う
            res = subprocess.run(cmd, capture_output=True, text=True, check=True)
            output = json.loads(res.stdout)

            if output.get("approved"):
                return True
            else:
                comments = output.get("comments", [])
                self._log_audit_failure(comments)
                return False
        except (subprocess.CalledProcessError, json.JSONDecodeError) as e:
            logger.error(f"Audit tool execution failed: {e}")
            # ツール自体の失敗は、とりあえずFalseとして扱うか、例外にするか
            return False
        except FileNotFoundError:
            logger.warning(
                f"Audit tool ({self.audit_tool_command[0]}) not found. Skipping audit."
            )
            # 本来はFailさせるべきだが、開発環境構築中なのでWarningに留めるケースも
            return False

    def _log_audit_failure(self, comments: list):
        with open(self.audit_log_path, "a") as f:
            f.write(f"\n## Audit Failed at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            for c in comments:
                f.write(f"- {c}\n")

    def run_uat_phase(self):
        """
        Phase 3.5: UATの生成と実行
        """
        if self.dry_run:
            logger.info("[DRY-RUN] Generating UAT with Playwright and running pytest...")
            return

        # 1. Generate UAT Code
        gen_cmd = self.ai_tool_command + [
            "--title",
            f"Generate UAT for Cycle {self.cycle_id}",
            "--description",
            f"Create Playwright tests in tests/e2e/ based on documents/CYCLE{self.cycle_id}/UAT.md",
        ]
        run_command(gen_cmd)

        # 2. Run Tests
        test_cmd = ["uv", "run", "pytest", "tests/e2e/"]
        try:
            run_command(test_cmd)
        except subprocess.CalledProcessError:
            self._trigger_fix("UAT Tests Failed. Please fix implementation or tests.")
            raise Exception("UAT Phase Failed")

    def finalize_cycle(self):
        """
        Phase 4: 自動マージ
        """
        if self.dry_run:
            logger.info("[DRY-RUN] Merging PR via gh CLI...")
            return

        cmd = ["gh", "pr", "merge", "--squash", "--delete-branch", "--admin"]
        run_command(cmd)
