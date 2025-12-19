from ac_cdd_core.config import settings


def test_config_agents_loaded():
    """Test that agents configuration is loaded correctly from files."""
    assert settings.agents.architect is not None
    # Check for content that actually exists in the prompts or defaults
    assert "Chief Systems Architect" in settings.agents.architect
    assert "Jules" in settings.agents.coder
    assert "Gemini" in settings.agents.auditor


def test_tools_config():
    """Test that tools configuration is loaded"""
    assert settings.tools.jules_cmd == "jules"
    assert settings.tools.gemini_cmd == "gemini"


def test_strict_audit_rules():
    """Verify stricter rules imply specific prompts"""
    # Verify based on actual content of tester.md
    assert "QA Engineer" in settings.agents.tester
