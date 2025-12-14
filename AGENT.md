# 🛠️ AGENT.md: Autonomous Development Template Refactoring Plan

## 1. Role & Objective (役割と目的)
あなたは、高度なPythonアプリケーション開発とAIエージェント統合のスペシャリストです。
現在、このリポジトリ (`autonomous-dev-template`) は、外部コマンドとして `jules` CLIツールを `subprocess` で呼び出す「CLIラッパー」として実装されています。
あなたのミッションは、このアーキテクチャを根本から見直し、**「Jules REST APIによるステートフルな対話」** と **「CLIによる互換性の維持」** を両立する **ハイブリッド・アーキテクチャ** へとリファクタリングすることです。

**Core Objective (最重要目標):**
環境変数 `JULES_API_KEY` の有無を検知し、自動的に以下の2つのモードを切り替える「Graceful Degradation（機能縮退）」メカニズムを実装してください。

1.  **API Mode (Priority)**: APIキーが存在する場合。REST APIを使用してセッションIDを維持し、継続的な対話と自律的な修正ループ（Human-in-the-loop）を実現する。
2.  **CLI Mode (Fallback)**: APIキーが存在しない場合。既存の `subprocess` 呼び出しを行い、ステートレスな単発実行機能を提供する（現状の機能を維持）。

---

## 2. Architecture & Design Pattern (アーキテクチャ設計)

この要件を満たすため、以下の設計パターンを採用して実装してください。

### 2.1. Abstract Base Class & Adapter Pattern
`scripts/ai_controller.py` 内のロジックが複雑化するのを防ぐため、エージェントとの通信部分を抽象化します。

* **`src/agent_interface.py` (新規作成)**:
    * 抽象基底クラス `AgentInterface` を定義。
    * 共通メソッド: `start_task(prompt)`, `send_message(prompt)`, `get_status()`.
* **`src/jules_api_client.py` (新規作成)**:
    * `AgentInterface` を継承。`httpx` を使用して Jules API (`https://jules.googleapis.com/v1alpha`) と通信。
    * **State Management**: セッションIDをローカルファイル（例: `.jules/session_state.json`）に永続化し、コンテキストを維持する責務を持つ。
* **`src/jules_cli_client.py` (新規作成)**:
    * `AgentInterface` を継承。既存の `subprocess.run` ロジックをカプセル化。
    * `send_message` は CLIではサポートされないため、内部で `start_task` (新規セッション作成) にフォールバックするか、警告を出して処理を継続する実装とする。

### 2.2. Factory Pattern
`scripts/ai_controller.py` の冒頭で、環境変数をチェックして適切なインスタンスを生成するファクトリ関数を実装してください。

```python
def get_agent_client() -> AgentInterface:
    api_key = os.getenv("JULES_API_KEY")
    if api_key:
        return JulesApiClient(api_key)
    else:
        typer.secho("⚠️ API Key not found. Running in CLI Mode (Stateless).", fg=typer.colors.YELLOW)
        return JulesCliClient()


ご提示いただいた「APIキーの有無によるGraceful Degradation（機能縮退）」の要件を組み込み、かつJulesが迷わず実装できるよう詳細な仕様、設計パターン、エッジケースの処理までを記述した、約4000文字規模の高密度プロンプトを作成しました。

この内容を AGENT.md として保存し、Julesに読み込ませてください。

AGENT.md （Julesへの指示書）
Markdown

# 🛠️ AGENT.md: Autonomous Development Template Refactoring Plan

## 1. Role & Objective (役割と目的)
あなたは、高度なPythonアプリケーション開発とAIエージェント統合のスペシャリストです。
現在、このリポジトリ (`autonomous-dev-template`) は、外部コマンドとして `jules` CLIツールを `subprocess` で呼び出す「CLIラッパー」として実装されています。
あなたのミッションは、このアーキテクチャを根本から見直し、**「Jules REST APIによるステートフルな対話」** と **「CLIによる互換性の維持」** を両立する **ハイブリッド・アーキテクチャ** へとリファクタリングすることです。

**Core Objective (最重要目標):**
環境変数 `JULES_API_KEY` の有無を検知し、自動的に以下の2つのモードを切り替える「Graceful Degradation（機能縮退）」メカニズムを実装してください。

1.  **API Mode (Priority)**: APIキーが存在する場合。REST APIを使用してセッションIDを維持し、継続的な対話と自律的な修正ループ（Human-in-the-loop）を実現する。
2.  **CLI Mode (Fallback)**: APIキーが存在しない場合。既存の `subprocess` 呼び出しを行い、ステートレスな単発実行機能を提供する（現状の機能を維持）。

---

## 2. Architecture & Design Pattern (アーキテクチャ設計)

この要件を満たすため、以下の設計パターンを採用して実装してください。

### 2.1. Abstract Base Class & Adapter Pattern
`scripts/ai_controller.py` 内のロジックが複雑化するのを防ぐため、エージェントとの通信部分を抽象化します。

* **`src/agent_interface.py` (新規作成)**:
    * 抽象基底クラス `AgentInterface` を定義。
    * 共通メソッド: `start_task(prompt)`, `send_message(prompt)`, `get_status()`.
* **`src/jules_api_client.py` (新規作成)**:
    * `AgentInterface` を継承。`httpx` を使用して Jules API (`https://jules.googleapis.com/v1alpha`) と通信。
    * **State Management**: セッションIDをローカルファイル（例: `.jules/session_state.json`）に永続化し、コンテキストを維持する責務を持つ。
* **`src/jules_cli_client.py` (新規作成)**:
    * `AgentInterface` を継承。既存の `subprocess.run` ロジックをカプセル化。
    * `send_message` は CLIではサポートされないため、内部で `start_task` (新規セッション作成) にフォールバックするか、警告を出して処理を継続する実装とする。

### 2.2. Factory Pattern
`scripts/ai_controller.py` の冒頭で、環境変数をチェックして適切なインスタンスを生成するファクトリ関数を実装してください。

```python
def get_agent_client() -> AgentInterface:
    api_key = os.getenv("JULES_API_KEY")
    if api_key:
        return JulesApiClient(api_key)
    else:
        typer.secho("⚠️ API Key not found. Running in CLI Mode (Stateless).", fg=typer.colors.YELLOW)
        return JulesCliClient()
3. Implementation Details (詳細実装要件)
A. Environment & Dependencies
Dependency Management:

pyproject.toml に httpx ライブラリを追加してください（uv add httpx コマンドを使用すること）。

Configuration:

.env ファイルの読み込みには既存の python-dotenv を利用します。

.jules/ ディレクトリ（セッション状態保存用）が .gitignore に含まれていることを確認し、なければ追加してください。

B. src/jules_api_client.py (The Brain)
このクラスはシステムの核となります。以下の仕様を厳守してください。

Session Persistence:

メソッド呼び出し時に、.jules/session_state.json から {"last_session_id": "..."} を読み込みます。

有効なセッションIDがあれば SendMessage エンドポイントを使用。

なければ CreateSession エンドポイントを使用し、返却されたIDを保存します。

Error Handling:

APIから 401 Unauthorized が返ってきた場合（キーが無効など）、例外を投げずに 動的にCLIモードへフォールバック するロジックを含めるとさらに理想的です（Optional）。

Logging:

API通信のログ（リクエスト/レスポンス）をデバッグモードで出力できるようにしてください。

C. Refactoring scripts/ai_controller.py
既存の typer コマンド群を、AgentInterface を使う形に書き換えます。

Command 1: strict_review (厳格レビュー)
共通フロー: git diff 取得 -> gemini CLIでレビュー -> JSON取得。

修正ロジックの分岐:

レビュー結果 (review_data) に has_issues: true が含まれる場合：

client.send_message(instructions) を呼び出す。

API Mode: 前回のセッション（コードを書いた本人）に対して「修正指示」が飛ぶため、文脈を理解した修正が行われる。

CLI Mode: 内部で jules remote new が走り、新規セッションとして修正が行われる。

Command 2: auto_fix (自動修正)
共通フロー: pytest 実行 -> 失敗ログ取得。

修正ロジックの分岐:

client.send_message(error_logs) を呼び出す。

API Mode: 「前回の修正は失敗しました。以下のログを見て再修正してください」という文脈が伝わる。これにより、Julesは「試行錯誤」が可能になる。

CLI Mode: 単にログを入力とした新規タスクが発行される。

Command 3: chat (新規追加)
機能: 開発者がターミナルから直接、現在のセッションに介入するためのコマンド。

Usage: uv run python scripts/ai_controller.py chat "もっと効率的なアルゴリズムにして"

実装:

API Mode: アクティブなセッションに SendMessage。

CLI Mode: 「CLIモードでは対話機能はサポートされていません」とエラー終了、または新規セッション作成にフォールバック。

4. Execution Plan (Julesへのタスクリスト)
以下の手順で実装を進めてください。各ステップで既存機能を破壊しないよう注意してください。

Setup: uv add httpx を実行し、pyproject.toml と uv.lock を更新する。

Skeleton: src/ 配下に agent_interface.py, jules_api_client.py, jules_cli_client.py の空ファイルを作成する。

Interface Definition: AgentInterface クラスを定義する。

CLI Adapter: 既存の ai_controller.py 内の subprocess ロジックを JulesCliClient に移植する。この時点で動作が変わらないことを確認する。

API Adapter: JulesApiClient を実装する。CreateSession, SendMessage のAPIコールとセッションID保存ロジックを実装する。

Controller Integration: scripts/ai_controller.py を修正し、Factory Patternを用いてクライアントを切り替えるようにする。

Logic Enhancement: strict_review, auto_fix コマンド内で、send_message を活用したループ処理ロジックを完成させる。

Verification: 最後に、.env のAPIキーを一時的に無効化（コメントアウト）しても、CLIモードとして動作することを確認するコードまたは手動テスト手順を追加する。

5. Constraints & Guidelines (制約事項)
Code Style: Python 3.12+ の機能を活用し、Type Hintingを必須とします。

Contracts: contracts/schemas.py は一切変更しないでください。

Safety: APIキーは絶対にログ出力や画面表示しないでください。

Robustness: APIモードでの通信エラー時は、ユーザーにわかりやすいエラーメッセージを表示し、可能であればCLIモードでの再試行を提案してください。

さあ、このリポジトリを「CLIツール」から「自律型AIパートナー」へと進化させましょう。実装を開始してください。


---

### 実行方法

この内容で `AGENT.md` を作成・コミットした後、以下のコマンドを実行してください。

```bash
jules remote new --repo <your-username>/autonomous-dev-template \
  --session "Architecture Upgrade: Hybrid Mode" \
  "リポジトリルートにある AGENT.md を読み、そこに記載されたハイブリッドアーキテクチャへのリファクタリングを実行してください。APIモードとCLIモードの自動切り替え機能の実装が最優先です。"
これで、JulesはAPIキーの有無を検知し、あなたの開発環境に合わせて賢く振る舞うようになります。