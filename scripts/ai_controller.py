import typer
import subprocess
import shutil
import sys
import json
import concurrent.futures
from pathlib import Path
from typing import Optional

app = typer.Typer(help="AI Controller for autonomous development workflows.")

def check_dependency(command: str, install_hint: str = ""):
    """Checks if a command exists in the PATH."""
    if not shutil.which(command):
        error_msg = f"Error: '{command}' command not found."
        if install_hint:
            error_msg += f" {install_hint}"
        typer.secho(error_msg, fg=typer.colors.RED)
        sys.exit(1)

def run_subprocess(cmd: list[str], check: bool = True, capture_output: bool = True) -> subprocess.CompletedProcess:
    """Runs a subprocess with error handling."""
    try:
        return subprocess.run(cmd, check=check, capture_output=capture_output, text=True)
    except subprocess.CalledProcessError as e:
        if check:
            typer.secho(f"Error executing command: {' '.join(cmd)}", fg=typer.colors.RED)
            if e.stderr:
                typer.secho(e.stderr, fg=typer.colors.RED)
            sys.exit(e.returncode)
        return e
    except Exception as e:
        typer.secho(f"Unexpected error: {str(e)}", fg=typer.colors.RED)
        sys.exit(1)

@app.command()
def strict_review():
    """
    Runs git diff, reviews it with Gemini, and triggers a Jules session for fixes.
    """
    check_dependency("git")
    check_dependency("gemini")
    check_dependency("jules")

    # Get git diff
    diff_process = run_subprocess(["git", "diff", "HEAD"], capture_output=True)
    diff_content = diff_process.stdout.strip()

    if not diff_content:
        typer.secho("No changes detected (git diff HEAD is empty).", fg=typer.colors.YELLOW)
        return

    prompt = (
        "あなたはGoogleのStaff Engineerです。入力されたdiffに対し、セキュリティ・パフォーマンス・可読性の観点から厳格なレビューを行い、"
        "具体的な修正指示のみをJulesへの命令形式で出力してください。\n\n"
        f"Diff:\n{diff_content}"
    )

    typer.secho("Requesting review from Gemini...", fg=typer.colors.BLUE)

    # Call Gemini
    gemini_process = run_subprocess(["gemini", "-p", prompt], capture_output=True)
    review_output = gemini_process.stdout.strip()

    if not review_output:
        typer.secho("Gemini returned no output.", fg=typer.colors.YELLOW)
        return

    typer.secho("Gemini review received. Starting Jules session...", fg=typer.colors.GREEN)

    # Start Jules session
    run_subprocess(["jules", "new", review_output], check=False, capture_output=False)


@app.command()
def auto_fix():
    """
    Runs tests, and if they fail, triggers a Jules session to fix them based on logs.
    """
    check_dependency("uv")
    check_dependency("jules")

    typer.secho("Running tests...", fg=typer.colors.BLUE)

    # Run pytest
    result = run_subprocess(["uv", "run", "pytest"], check=False, capture_output=True)

    if result.returncode == 0:
        typer.secho("All tests passed.", fg=typer.colors.GREEN)
        return

    typer.secho("Tests failed. Analyzing logs...", fg=typer.colors.RED)

    # Prepare prompt with logs
    logs = result.stdout + "\n" + result.stderr
    prompt = f"テストが失敗しました。ログを分析し修正してください。\n\nLogs:\n{logs}"

    typer.secho("Starting Jules session for auto-fix...", fg=typer.colors.BLUE)
    run_subprocess(["jules", "new", prompt], check=False, capture_output=False)


@app.command()
def triage_issues():
    """
    Fetches assigned GitHub issues and starts a Jules session for each.
    """
    check_dependency("gh", "Please install GitHub CLI.")
    check_dependency("jules")

    typer.secho("Fetching assigned issues...", fg=typer.colors.BLUE)

    # Fetch issues
    cmd = ["gh", "issue", "list", "--assignee", "@me", "--json", "number,title,body"]
    result = run_subprocess(cmd, capture_output=True)

    try:
        issues = json.loads(result.stdout)
    except json.JSONDecodeError:
        typer.secho("Failed to parse GitHub issues JSON.", fg=typer.colors.RED)
        sys.exit(1)

    if not issues:
        typer.secho("No assigned issues found.", fg=typer.colors.YELLOW)
        return

    typer.secho(f"Found {len(issues)} issues. Starting Jules sessions...", fg=typer.colors.GREEN)

    for issue in issues:
        number = issue.get("number")
        title = issue.get("title")
        body = issue.get("body")

        prompt = f"Issue #{number}: {title}\n\n{body}"
        typer.secho(f"Processing Issue #{number}...", fg=typer.colors.BLUE)

        run_subprocess(["jules", "new", prompt], check=False, capture_output=False)


@app.command()
def implement_plan(
    file: Path = typer.Argument(
        Path("documents/TODO.md"),
        help="Path to the task list file.",
        exists=True,
        dir_okay=False,
        readable=True
    ),
    parallel: int = typer.Option(5, help="Number of parallel tasks.")
):
    """
    Reads a task list file and executes Jules sessions in parallel for each line.
    """
    check_dependency("jules")

    tasks = [line.strip() for line in file.read_text().splitlines() if line.strip()]

    if not tasks:
        typer.secho("No tasks found in file.", fg=typer.colors.YELLOW)
        return

    typer.secho(f"Found {len(tasks)} tasks. Executing with {parallel} threads...", fg=typer.colors.GREEN)

    def execute_task(task_prompt):
        typer.secho(f"Starting task: {task_prompt[:50]}...", fg=typer.colors.BLUE)
        try:
            subprocess.run(["jules", "new", task_prompt], check=True, capture_output=False)
            return f"Completed: {task_prompt[:30]}..."
        except subprocess.CalledProcessError as e:
            return f"Failed: {task_prompt[:30]}..."

    with concurrent.futures.ThreadPoolExecutor(max_workers=parallel) as executor:
        future_to_task = {executor.submit(execute_task, task): task for task in tasks}
        for future in concurrent.futures.as_completed(future_to_task):
            task = future_to_task[future]
            try:
                data = future.result()
                typer.secho(data, fg=typer.colors.GREEN)
            except Exception as exc:
                typer.secho(f"Task generated an exception: {exc}", fg=typer.colors.RED)

if __name__ == "__main__":
    app()
