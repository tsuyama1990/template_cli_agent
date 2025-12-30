from unittest.mock import AsyncMock, MagicMock, patch

import pytest
from ac_cdd_core.cli import app
from ac_cdd_core.state import CycleState
from typer.testing import CliRunner


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def mock_project_manager():
    with patch("ac_cdd_core.services.project.ProjectManager") as mock:
        mock_instance = MagicMock()
        mock.return_value = mock_instance
        yield mock_instance


@pytest.fixture
def mock_session_manager():
    with patch("ac_cdd_core.session_manager.SessionManager") as mock:
        mock.SESSION_FILE = MagicMock()
        mock.SESSION_FILE.exists.return_value = True
        mock.load_session.return_value = {
            "session_id": "test-session",
            "integration_branch": "dev/test-session",
        }
        mock.get_integration_branch.return_value = "dev/test-session"
        # Mock async method - actually load_or_reconcile_session is now SYNC in source code
        mock.load_or_reconcile_session.return_value = ("test-session", "dev/test-session")
        yield mock


@pytest.fixture
def mock_session_validator():
    with patch("ac_cdd_core.cli.SessionValidator") as mock:
        mock_instance = AsyncMock()
        mock.return_value = mock_instance
        mock_instance.validate.return_value = (True, "")
        yield mock_instance


@pytest.fixture
def mock_graph_builder():
    with patch("ac_cdd_core.cli.GraphBuilder") as mock:
        mock_instance = MagicMock()  # GraphBuilder itself is sync
        mock.return_value = mock_instance
        # Mock ainvoke to return a state dict
        # The graph compiled object has ainvoke
        mock_graph = MagicMock()
        mock_instance.build_architect_graph.return_value = mock_graph
        mock_graph.ainvoke = AsyncMock(
            return_value=CycleState(
                session_id="test-session", integration_branch="dev/test-session"
            )
        )

        mock_coder_graph = MagicMock()
        mock_instance.build_coder_graph.return_value = mock_coder_graph
        mock_coder_graph.ainvoke = AsyncMock(return_value=CycleState())
        yield mock_instance


@pytest.fixture
def mock_git_manager():
    # Since imports are inside functions, we patch ac_cdd_core.services.git_ops.GitManager
    with patch("ac_cdd_core.services.git_ops.GitManager") as mock:
        mock_instance = AsyncMock()
        mock.return_value = mock_instance
        mock_instance.create_final_pr.return_value = "https://github.com/user/repo/pull/1"
        yield mock_instance


def test_init_command(runner, mock_project_manager):
    """Test init command execution."""
    # Mock check_environment to pass
    with patch("ac_cdd_core.cli.check_environment"):
        result = runner.invoke(app, ["init"])
        assert result.exit_code == 0
        assert "Initialization Complete" in result.stdout
        mock_project_manager.initialize_project.assert_called_once()


def test_gen_cycles_command(runner, mock_graph_builder, mock_session_manager):
    """Test gen-cycles command execution."""
    with patch("ac_cdd_core.messages.ensure_api_key"):
        result = runner.invoke(app, ["gen-cycles", "--session", "test-session"])
        assert result.exit_code == 0
        assert "Architect Phase Complete" in result.stdout
        mock_session_manager.save_session.assert_called()


def test_run_cycle_command(
    runner, mock_graph_builder, mock_session_manager, mock_session_validator
):
    """Test run-cycle command execution."""
    with patch("ac_cdd_core.messages.ensure_api_key"):
        result = runner.invoke(app, ["run-cycle", "--id", "01"])
        assert result.exit_code == 0
        assert "Cycle 01 Implementation Request Sent" in result.stdout


def test_finalize_session_command(runner, mock_session_manager, mock_git_manager):
    """Test finalize-session command."""
    with patch(
        "ac_cdd_core.session_manager.SessionManager.validate_session", return_value=(True, "")
    ):
        result = runner.invoke(app, ["finalize-session", "--session", "test-session"])
        assert result.exit_code == 0
        assert "Final PR Created" in result.stdout
        mock_git_manager.create_final_pr.assert_called()


def test_check_environment_missing_keys():
    """Test check_environment detects missing API keys."""
    with (
        patch("ac_cdd_core.utils.check_api_key", return_value=False),
        patch("rich.console.Console.print") as _mock_print,
    ):
        from ac_cdd_core.cli import check_environment

        with pytest.raises(SystemExit):
            check_environment()


def test_check_environment_all_present():
    """Test check_environment passes when all keys present."""
    with (
        patch("ac_cdd_core.utils.check_api_key", return_value=True),
        patch("subprocess.run", return_value=MagicMock(returncode=0)),
    ):
        from ac_cdd_core.cli import check_environment

        # Should not raise
        check_environment()
