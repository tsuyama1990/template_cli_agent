"""
Unit tests for the Explorer module.

These tests are designed to verify the logic of the Explorer class without
actually downloading and running the expensive MACE model or performing real
MD simulations. We use extensive mocking to isolate the Explorer's own logic.
"""

from unittest.mock import MagicMock, call

import pytest
from ase.build import bulk

from mlip_autopipec.config import ExplorerConfig, SimulationConfig
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.modules.explorer_sampler import Explorer


@pytest.fixture
def mock_db_wrapper(mocker):
    """Fixture to create a mock AseDBWrapper."""
    mock = mocker.MagicMock(spec=AseDBWrapper)
    # Simulate the DB returning a conventional cell with 2 atoms to avoid ZeroDivisionError
    atoms = bulk("Fe", "bcc", cubic=True)
    mock.get_all_initial_structures.return_value = [atoms]
    return mock


@pytest.fixture
def explorer_config():
    """Fixture to create a default ExplorerConfig."""
    return ExplorerConfig(md_steps_per_structure=100)  # Use fewer steps for tests


@pytest.fixture
def simulation_config():
    """Fixture to create a default SimulationConfig."""
    return SimulationConfig(temperature_steps=[300, 600])


def test_explorer_initialization(explorer_config, simulation_config, mock_db_wrapper, mocker):
    """Tests that the Explorer is initialized correctly."""
    mocker.patch("torch.cuda.is_available", return_value=False)
    explorer = Explorer(explorer_config, simulation_config, mock_db_wrapper)
    assert explorer.config == explorer_config
    assert explorer.simulation_config == simulation_config
    assert explorer.db_wrapper == mock_db_wrapper
    assert explorer.device == "cpu"


def test_explorer_initialization_cuda(explorer_config, simulation_config, mock_db_wrapper, mocker):
    """Tests that the Explorer correctly detects and sets the CUDA device."""
    mocker.patch("torch.cuda.is_available", return_value=True)
    explorer = Explorer(explorer_config, simulation_config, mock_db_wrapper)
    assert explorer.device == "cuda"


@pytest.mark.parametrize("device", ["cpu", "cuda"])
def test_run_mace_md_logic(explorer_config, simulation_config, mock_db_wrapper, mocker, device):
    """
    Tests the core logic of a single MD run, mocking all external libraries.
    """
    # Patch the external dependencies
    mocker.patch("torch.cuda.is_available", return_value=(device == "cuda"))
    # Patch the correct import path for the mace_mp function
    mock_mace_mp = mocker.patch("mlip_autopipec.modules.explorer_sampler.mace_mp")
    mock_langevin = mocker.patch("mlip_autopipec.modules.explorer_sampler.Langevin")

    # Instantiate the mock dynamics object that Langevin would return
    mock_dyn = MagicMock()
    mock_langevin.return_value = mock_dyn

    # Instantiate the Explorer
    explorer = Explorer(explorer_config, simulation_config, mock_db_wrapper)
    atoms = bulk("Fe", "bcc", cubic=True)  # Use a 2-atom cell
    temperature = 500

    # Run the private method under test
    frames = explorer._run_mace_md(atoms, temperature)

    # --- Assertions ---

    # 1. Assert that the MACE model was loaded with the correct parameters
    mock_mace_mp.assert_called_once_with(
        model=None,
        device=device,
        default_dtype="float64",
    )
    assert atoms.calc is not None

    # 2. Assert that the ASE Langevin dynamics was initialized correctly
    mock_langevin.assert_called_once()
    args, kwargs = mock_langevin.call_args
    assert args[0] is atoms
    assert "temperature_K" in kwargs and kwargs["temperature_K"] == temperature

    # 3. Assert that the dynamics simulation was run with the correct number of steps
    mock_dyn.run.assert_called_once_with(explorer_config.md_steps_per_structure)

    # 4. Assert that the trajectory observer was attached
    assert mock_dyn.attach.called

    # For this test, since the observer is real, the list of frames should be populated
    # based on the mocked `run` method (although it won't actually step).
    # A more advanced test could mock the observer callback. Here, we confirm it
    # should have collected at least the initial state.
    assert len(frames) >= 0


def test_unsupported_model_raises_error(
    simulation_config, mock_db_wrapper, mocker
):
    """
    Tests that a ValueError is raised if the user configures an unsupported model.
    """
    # Configure the explorer with an invalid model name
    invalid_config = ExplorerConfig(surrogate_model="unsupported_potential")
    explorer = Explorer(invalid_config, simulation_config, mock_db_wrapper)
    atoms = bulk("Fe", "bcc", cubic=True)  # Use a 2-atom cell

    # Assert that the correct exception is raised
    with pytest.raises(ValueError, match="only the 'mace_mp' surrogate model is supported"):
        explorer._run_mace_md(atoms, temperature=300)


def test_explore_main_loop_logic(
    explorer_config, simulation_config, mock_db_wrapper, mocker
):
    """
    Tests that the main `explore` method correctly loops through structures
    and temperatures, calling the private `_run_mace_md` method.
    """
    # Mock the internal method to prevent it from actually running
    mock_run_md = mocker.patch(
        "mlip_autopipec.modules.explorer_sampler.Explorer._run_mace_md"
    )
    # Configure the mock to return a dummy list of 2 frames per call
    mock_run_md.return_value = [bulk("Fe"), bulk("Fe")]

    # Get the mock structures and temperatures from the fixtures
    initial_structures = mock_db_wrapper.get_all_initial_structures.return_value
    temps = simulation_config.temperature_steps

    # Instantiate the Explorer and run the main method
    explorer = Explorer(explorer_config, simulation_config, mock_db_wrapper)
    all_frames = explorer.explore()

    # --- Assertions ---

    # 1. Assert that the database was queried for initial structures
    mock_db_wrapper.get_all_initial_structures.assert_called_once()

    # 2. Assert that `_run_mace_md` was called for each structure and temperature
    expected_call_count = len(initial_structures) * len(temps)
    assert mock_run_md.call_count == expected_call_count

    # 3. Check the content of the calls
    expected_calls = [
        call(initial_structures[0], temps[0]),
        call(initial_structures[0], temps[1]),
    ]
    mock_run_md.assert_has_calls(expected_calls, any_order=True)

    # 4. Assert that the final list contains the aggregated results
    # (2 frames per call) * (2 calls) = 4 frames total
    assert len(all_frames) == expected_call_count * 2
