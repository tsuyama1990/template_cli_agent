import pytest
from pydantic import ValidationError

from mlip_autopipec.config import ConfigExpander, UserInputConfig


@pytest.fixture
def fept_input() -> UserInputConfig:
    """Returns a user input for a FePt alloy."""
    return UserInputConfig.model_validate(
        {
            "system": {"elements": ["Fe", "Pt"], "composition": "FePt"},
        }
    )


@pytest.fixture
def si_input() -> UserInputConfig:
    """Returns a user input for Silicon."""
    return UserInputConfig.model_validate(
        {
            "system": {"elements": ["Si"], "composition": "Si"},
        }
    )


@pytest.fixture
def temp_input() -> UserInputConfig:
    """Returns a user input with a specified temperature range."""
    return UserInputConfig.model_validate(
        {
            "system": {"elements": ["Al"], "composition": "Al"},
            "simulation": {"temperature": [300, 1000]},
        }
    )


def test_expand_magnetic_alloy(fept_input: UserInputConfig):
    """Tests the expansion of a magnetic FePt alloy user input."""
    expander = ConfigExpander()
    full_config = expander.expand(fept_input)

    # System assertions
    assert full_config.system.structure_type == "alloy"
    assert full_config.system.melting_point_guess > 0

    # DFT assertions
    assert full_config.dft_compute.magnetism == "ferromagnetic"
    assert full_config.dft_compute.ecutwfc == 90.0
    assert full_config.dft_compute.ecutrho == 720.0


def test_expand_nonmagnetic_covalent(si_input: UserInputConfig):
    """Tests the expansion of a non-magnetic Si covalent user input."""
    expander = ConfigExpander()
    full_config = expander.expand(si_input)

    # System assertions
    assert full_config.system.structure_type == "covalent"

    # DFT assertions
    assert full_config.dft_compute.magnetism is None
    assert full_config.dft_compute.ecutwfc == 40.0
    assert full_config.dft_compute.ecutrho == 320.0


def test_temperature_interpolation(temp_input: UserInputConfig):
    """Tests that a user-provided temperature range is interpolated."""
    expander = ConfigExpander()
    full_config = expander.expand(temp_input)

    # Simulation assertions
    assert full_config.simulation.temperature_steps == [300, 650, 1000]


def test_default_temperature_generation(fept_input: UserInputConfig):
    """Tests that a default temperature range is generated if none is provided."""
    expander = ConfigExpander()
    full_config = expander.expand(fept_input)

    # Simulation assertions
    assert len(full_config.simulation.temperature_steps) > 1
    assert full_config.simulation.temperature_steps[0] == 300
    assert (
        full_config.simulation.temperature_steps[-1] <= 0.8 * full_config.system.melting_point_guess
    )


def test_invalid_input_raises_error():
    """Tests that invalid user input raises a ValidationError."""
    with pytest.raises(ValidationError):
        UserInputConfig.model_validate({"system": {"elements": []}})  # No elements


def test_missing_pseudopotential_raises_error(fept_input: UserInputConfig):
    """
    Tests that a FullConfig validation error is raised if an element
    is missing a pseudopotential.
    """
    expander = ConfigExpander()
    full_config = expander.expand(fept_input)

    # Manually create an invalid state
    full_config.dft_compute.pseudopotentials.pop("Pt")

    with pytest.raises(ValidationError) as excinfo:
        # Re-validate the modified model to trigger the validator
        full_config.model_validate(full_config.model_dump())

    assert "Missing pseudopotentials for elements: ['Pt']" in str(excinfo.value)
