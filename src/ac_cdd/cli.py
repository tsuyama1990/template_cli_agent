import asyncio
import os
import sys

import logfire
import typer
from rich.console import Console

from .config import settings
from .orchestrator import CycleOrchestrator
from .tools import ToolWrapper, ToolNotFoundError
from .agents import AgentDeps, auditor_agent, coder_agent
from .domain_models import AuditResult, FileArtifact
from .utils import logger

# Initialize Logfire
# We should only configure if we want it to send data, usually needs a token.
# If no token, it might warn or default to console.
# Ideally we load this conditionally, but the requirement was "import logfire; logfire.configure()"
# We suppress the error if not authenticated to allow local usage without logfire auth
try:
    logfire.configure()
except Exception as e:
    # Fallback or silent ignore if strict requirements say "ensure observability" but don't want to crash.
    # The error is LogfireConfigError.
    # We can try to configure with send_to_logfire=False to just get console logs if that's the intent
    # or just suppress it.
    # Let's try to configure for console only or just pass.
    # logfire.configure(send_to_logfire='if-token-present') seems better if available,
    # but based on docs or common pattern:
    import logging
    logging.getLogger("ac_cdd").warning(f"Logfire not configured: {e}")

app = typer.Typer(
    name="AC-CDD CLI",
    help="Autonomous Cycle-Based Contract-Driven Development Environment",
    add_completion=False,
)
console = Console()


@app.command()
def new_cycle(cycle_id: str = typer.Argument(..., help="Cycle ID (e.g. '01')")):
    """
    Creates a new cycle directory and scaffolds files.
    """
    import shutil
    from pathlib import Path

    base_dir = Path(settings.paths.documents_dir)
    cycle_dir = base_dir / f"CYCLE{cycle_id}"

    if cycle_dir.exists():
        console.print(f"[bold red]Error:[/bold red] {cycle_dir} already exists.")
        raise typer.Exit(code=1)

    cycle_dir.mkdir(parents=True)
    console.print(f"[green]Created {cycle_dir}[/green]")

    # Copy templates if they exist
    templates_dir = Path(settings.paths.templates) / "cycle"
    if templates_dir.exists():
        for item in ["SPEC.md", "UAT.md", "schema.py"]:
            src = templates_dir / item
            dst = cycle_dir / item
            if src.exists():
                shutil.copy(src, dst)
                console.print(f"  - Copied {item}")
    else:
        # Create empty files if templates don't exist
        for item in ["SPEC.md", "UAT.md", "schema.py"]:
            (cycle_dir / item).touch()
            console.print(f"  - Created empty {item}")


@app.command()
def start_cycle(
    cycle_ids: list[str] = typer.Argument(..., help="List of Cycle IDs to execute"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Simulate execution"),
    auto_next: bool = typer.Option(False, "--auto-next", help="Automatically prepare next cycle"),
):
    """
    Starts the AC-CDD automation cycle for the given ID(s).
    """
    async def _run():
        for cycle_id in cycle_ids:
            console.print(f"[bold blue]=== Starting Cycle {cycle_id} ===[/bold blue]")
            orchestrator = CycleOrchestrator(cycle_id, dry_run=dry_run, auto_next=auto_next)

            # Simple progress reporting could be added here
            await orchestrator.execute_all()

            console.print(f"[bold green]=== Cycle {cycle_id} Completed ===[/bold green]")

    asyncio.run(_run())


@app.command()
def audit(
    target: str = typer.Option("HEAD", "--target", help="Git ref or file path to audit"),
    fix: bool = typer.Option(False, "--fix", help="Auto-apply fixes"),
):
    """
    Audits the current changes using Gemini (via Auditor Agent).
    """
    async def _run_audit():
        from pathlib import Path
        from .tools import GitClient

        console.print(f"[bold]Auditing target: {target}[/bold]")

        # 1. Get Diff or Content
        content_to_audit = ""
        if Path(target).exists():
            content_to_audit = Path(target).read_text()
        else:
            # Assume git ref
            try:
                git = GitClient()
                content_to_audit = git.get_diff(target)
            except Exception as e:
                console.print(f"[red]Failed to get diff: {e}[/red]")
                raise typer.Exit(1)

        if not content_to_audit.strip():
            console.print("[yellow]No content to audit.[/yellow]")
            return

        # 2. Call Auditor Agent
        console.print("[cyan]Sending to Auditor Agent...[/cyan]")
        deps = AgentDeps(documents_dir=Path(settings.paths.documents_dir))

        # We wrap content in a prompt
        user_task = f"Audit the following code/diff:\n\n{content_to_audit}"

        try:
            result = await auditor_agent.run(user_task, deps=deps)
            audit_res: AuditResult = result.data

            if audit_res.is_approved:
                console.print("[bold green]Audit Passed![/bold green]")
            else:
                console.print("[bold red]Audit Failed[/bold red]")
                for issue in audit_res.critical_issues:
                    console.print(f"- [Red]Issue[/Red]: {issue}")
                for sug in audit_res.suggestions:
                    console.print(f"- [Blue]Suggestion[/Blue]: {sug}")

                if fix:
                    console.print("[bold blue]Applying Fixes...[/bold blue]")
                    fix_instructions = (
                        f"Fix the following issues in the code:\n" +
                        "\n".join(audit_res.critical_issues + audit_res.suggestions) +
                        f"\n\nOriginal Code/Diff:\n{content_to_audit}"
                    )
                    fix_result = await coder_agent.run(fix_instructions, deps=deps)
                    files: list[FileArtifact] = fix_result.data
                    for f in files:
                        p = Path(f.path)
                        p.parent.mkdir(parents=True, exist_ok=True)
                        p.write_text(f.content)
                        console.print(f"  - Updated {p}")

        except Exception as e:
            console.print(f"[bold red]Error during audit: {e}[/bold red]")
            raise typer.Exit(1)

    asyncio.run(_run_audit())


@app.command()
def fix(
    prompt: str = typer.Argument(..., help="Instructions for the fix"),
    file: str = typer.Option(None, "--file", help="Specific file context"),
):
    """
    Applies a fix using the Coder Agent.
    """
    async def _run_fix():
        from pathlib import Path

        context = ""
        if file and Path(file).exists():
            context = f"\n\nTarget File ({file}):\n{Path(file).read_text()}"

        full_prompt = prompt + context
        deps = AgentDeps(documents_dir=Path(settings.paths.documents_dir))

        console.print("[cyan]Asking Coder Agent to fix...[/cyan]")
        try:
            result = await coder_agent.run(full_prompt, deps=deps)
            files: list[FileArtifact] = result.data

            for f in files:
                p = Path(f.path)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(f.content)
                console.print(f"[green]Applied fix to {p}[/green]")

        except Exception as e:
            console.print(f"[red]Fix failed: {e}[/red]")
            raise typer.Exit(1)

    asyncio.run(_run_fix())


@app.command()
def doctor():
    """
    Checks the environment for required tools.
    """
    active_tools = [
        settings.tools.gh_cmd,
        settings.tools.audit_cmd,
        settings.tools.uv_cmd,
        settings.tools.mypy_cmd,
        "ruff"
    ]

    console.print("[bold]Checking Environment...[/bold]")
    all_ok = True
    import shutil

    for tool in active_tools:
        if shutil.which(tool):
            console.print(f"[green]✔ {tool} found[/green]")
        else:
            console.print(f"[red]✘ {tool} NOT found[/red]")
            all_ok = False

    # Check API Key
    if "GEMINI_API_KEY" not in os.environ:
        console.print("[red]✘ GEMINI_API_KEY not set[/red]")
        all_ok = False
    else:
        console.print("[green]✔ GEMINI_API_KEY set[/green]")

    if all_ok:
        console.print("[bold green]System is healthy.[/bold green]")
    else:
        console.print("[bold red]System has issues. Please install missing tools.[/bold red]")
        raise typer.Exit(1)


def main():
    app()

if __name__ == "__main__":
    main()
