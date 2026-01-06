from collections.abc import Iterator
from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from ac_cdd_core.cli import app
from ac_cdd_core.domain_models import CycleManifest, ProjectManifest
from typer.testing import CliRunner

runner = CliRunner()


@pytest.fixture
def mock_dependencies() -> Iterator[None]:
    with (
        patch("ac_cdd_core.cli.utils.check_api_key", return_value=True),
        patch("shutil.which", return_value="/usr/bin/git"),
        patch("ac_cdd_core.cli.ProjectManager"),
        patch("ac_cdd_core.cli.SessionManager"),
        patch("ac_cdd_core.cli._WorkflowServiceHolder.get"),
    ):
        yield


def test_init_command(mock_dependencies: None) -> None:
    # init calls SessionManager().load_manifest() and possibly git.ensure_state_branch
    with patch("ac_cdd_core.cli.SessionManager") as mock_sm_cls:
        mock_mgr = mock_sm_cls.return_value
        # load_manifest is async, but init runs asyncio.run(_init_state())
        mock_mgr.load_manifest = AsyncMock(return_value=None)
        mock_mgr.git.ensure_state_branch = AsyncMock()

        result = runner.invoke(app, ["init"])
        assert result.exit_code == 0
        assert "Initialization Complete" in result.stdout


def test_gen_cycles_command(mock_dependencies: None) -> None:
    mock_workflow = MagicMock()
    mock_workflow.run_gen_cycles = AsyncMock()
    with patch("ac_cdd_core.cli._WorkflowServiceHolder.get", return_value=mock_workflow):
        result = runner.invoke(app, ["gen-cycles", "--cycles", "3"])
        assert result.exit_code == 0
        mock_workflow.run_gen_cycles.assert_awaited_once_with(3, None)


def test_gen_cycles_command_handles_exception(mock_dependencies: None) -> None:
    """Tests that the gen-cycles command handles exceptions gracefully."""
    error_message = "A critical generation error occurred."
    with patch("asyncio.run", side_effect=Exception(error_message)):
        result = runner.invoke(app, ["gen-cycles", "--cycles", "3"])

    assert result.exit_code == 1
    assert "An unexpected error occurred" in result.stdout
    assert error_message in result.stdout
    assert "Traceback" not in result.stdout


def test_run_cycle_command(mock_dependencies: None) -> None:
    mock_workflow = MagicMock()
    mock_workflow.run_cycle = AsyncMock()
    with patch("ac_cdd_core.cli._WorkflowServiceHolder.get", return_value=mock_workflow):
        result = runner.invoke(app, ["run-cycle", "--id", "01"])
        assert result.exit_code == 0
        mock_workflow.run_cycle.assert_awaited_once_with(
            cycle_id="01", resume=False, auto=True, start_iter=1, project_session_id=None
        )


def test_run_cycle_command_handles_exception(mock_dependencies: None) -> None:
    """Tests that the CLI's global exception handler catches errors gracefully."""
    error_message = "A critical workflow error occurred."
    with patch("asyncio.run", side_effect=Exception(error_message)):
        result = runner.invoke(app, ["run-cycle", "--id", "01"])

    assert result.exit_code == 1
    assert "An unexpected error occurred" in result.stdout
    assert error_message in result.stdout
    assert "Traceback" not in result.stdout


def test_finalize_session_command(mock_dependencies: None) -> None:
    mock_workflow = MagicMock()
    mock_workflow.finalize_session = AsyncMock()
    with patch("ac_cdd_core.cli._WorkflowServiceHolder.get", return_value=mock_workflow):
        result = runner.invoke(app, ["finalize-session"])
        assert result.exit_code == 0
        mock_workflow.finalize_session.assert_awaited_once_with(None)


def test_list_actions_no_session(mock_dependencies: None) -> None:
    with patch("ac_cdd_core.cli.SessionManager") as mock_sm_cls:
        mock_mgr = mock_sm_cls.return_value
        mock_mgr.load_manifest = AsyncMock(return_value=None)

        result = runner.invoke(app, ["list-actions"])
        assert result.exit_code == 0
        assert "No active session found" in result.stdout


def test_list_actions_active_session(mock_dependencies: None) -> None:
    with patch("ac_cdd_core.cli.SessionManager") as mock_sm_cls:
        mock_mgr = mock_sm_cls.return_value
        manifest = ProjectManifest(
            project_session_id="p1",
            integration_branch="dev/p1",
            cycles=[CycleManifest(id="01", status="planned")],
        )
        mock_mgr.load_manifest = AsyncMock(return_value=manifest)

        result = runner.invoke(app, ["list-actions"])
        assert result.exit_code == 0
        assert "Active Session: p1" in result.stdout
        assert "run-cycle --id 01" in result.stdout
