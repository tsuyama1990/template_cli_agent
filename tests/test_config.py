from ac_cdd_core.config import settings


def test_config_agents_loaded():
    """Test that agents configuration is loaded correctly from files."""
    # Updated to reflect actual AgentsConfig fields: auditor, qa_analyst
    assert settings.agents.auditor is not None
    assert (
        "DEFAULT_AUDITOR_PROMPT" in settings.agents.auditor or "Auditor" in settings.agents.auditor
    )
    assert settings.agents.qa_analyst is not None
    assert (
        "DEFAULT_QA_ANALYST_PROMPT" in settings.agents.qa_analyst
        or "QA" in settings.agents.qa_analyst
    )


def test_tools_config():
    """Test that tools configuration is loaded"""
    assert settings.tools.jules_cmd == "jules"
    assert settings.tools.gemini_cmd == "gemini"
