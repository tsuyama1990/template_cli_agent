from unittest.mock import MagicMock, patch

import pytest
from ac_cdd_core.cli import app
from typer.testing import CliRunner

runner = CliRunner()


@pytest.fixture
def mock_dependencies():
    with (
        patch("ac_cdd_core.cli.check_environment"),
        patch("ac_cdd_core.services.project.ProjectManager"),
        patch("ac_cdd_core.graph.GraphBuilder") as mock_graph_cls,
        patch("ac_cdd_core.session_manager.SessionManager") as mock_session_mgr,
        patch("ac_cdd_core.messages.ensure_api_key"),
    ):
        mock_graph_instance = mock_graph_cls.return_value
        mock_graph_instance.build_architect_graph.return_value = MagicMock()
        mock_graph_instance.build_coder_graph.return_value = MagicMock()
        mock_graph_instance.cleanup = MagicMock()  # async?

        yield {"graph": mock_graph_instance, "session": mock_session_mgr}


def test_init_command(mock_dependencies):
    result = runner.invoke(app, ["init"])
    assert result.exit_code == 0
    assert "Initialization Complete" in result.stdout


def test_gen_cycles_command(mock_dependencies):
    # Mock async graph run
    mock_graph = mock_dependencies["graph"]
    # ainvoke returns final state
    mock_graph.build_architect_graph.return_value.ainvoke = MagicMock(
        side_effect=lambda state, *args, **kwargs: MagicMock(
            session_id="test-session",
            integration_branch="dev/test",
            get=lambda k: None,  # error is None
        )
    )

    # NOTE: typer using asyncio.run under the hood?
    # ac_cdd_core.cli.gen_cycles uses asyncio.run(_run())
    # So our mocks (which are Sync MagicMocks) might cause issues if awaited?
    # Python 3.8+ AsyncMock is needed for awaitables.
    from unittest.mock import AsyncMock

    mock_graph.build_architect_graph.return_value.ainvoke = AsyncMock(
        return_value=MagicMock(
            session_id="test-session", integration_branch="dev/test", get=lambda k: None
        )
    )
    mock_graph.cleanup = AsyncMock()

    result = runner.invoke(app, ["gen-cycles"])
    assert result.exit_code == 0
    assert "Architect Phase" in result.stdout
    assert "Using saved session" not in result.stdout  # Should be new


def test_run_cycle_command(mock_dependencies):
    mock_graph = mock_dependencies["graph"]
    mock_session = mock_dependencies["session"]
    from unittest.mock import AsyncMock

    mock_session.load_or_reconcile_session.return_value = {
        "session_id": "test",
        "integration_branch": "dev/test",
    }
    mock_session.validate_session.return_value = (True, None)

    mock_graph.build_coder_graph.return_value.ainvoke = AsyncMock(
        return_value={"current_phase": "complete"}
    )
    mock_graph.cleanup = AsyncMock()  # Fix await

    # Patch SessionValidator because run_cycle instantiates it
    with patch("ac_cdd_core.validators.SessionValidator") as mock_val_cls:
        mock_val_cls.return_value.validate = AsyncMock(return_value=(True, None))

        result = runner.invoke(app, ["run-cycle", "--id", "01"])

        if result.exit_code != 0:
            print(result.stdout)

        assert result.exit_code == 0
        assert "Cycle 01 Implementation Request Sent" in result.stdout


def test_finalize_session_command(mock_dependencies):
    mock_session = mock_dependencies["session"]
    mock_session.load_or_reconcile_session.return_value = {
        "session_id": "test",
        "integration_branch": "dev/test",
    }
    mock_session.validate_session.return_value = (True, None)

    with patch("ac_cdd_core.services.git_ops.GitManager") as mock_git_cls:
        mock_git_instance = mock_git_cls.return_value
        from unittest.mock import AsyncMock

        mock_git_instance.create_final_pr = AsyncMock(return_value="https://pr")

        result = runner.invoke(app, ["finalize-session"])

        assert result.exit_code == 0
        assert "Session Finalized" in result.stdout
