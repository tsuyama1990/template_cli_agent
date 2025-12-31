"""Unit tests for the LabellingService."""

import pytest
from unittest.mock import MagicMock
from ase import Atoms
import numpy as np

from mlip_autopipec.application.services import LabellingService
from mlip_autopipec.domain.models import DFTResult
from mlip_autopipec.config import Settings

@pytest.fixture
def mock_db_port():
    return MagicMock()

@pytest.fixture
def mock_runner_port():
    return MagicMock()

@pytest.fixture
def mock_settings():
    return Settings(dft_command="test cmd", pseudo_dir=".", ecutwfc=50.0)

SAMPLE_SUCCESS_OUTPUT = """
!    total energy              =     -11.45567302 Ry

Forces acting on atoms (cartesian axes, Ry/au):

     atom    1   force =    -0.00000101    0.00000000   -0.00000000

total stress  (Ry/bohr**3)                (kbar)     P=      -0.03
 0.00000004   -0.00000000    0.00000000
 0.00000000    0.00000004    0.00000000
 0.00000000    0.00000000    0.00000004
"""

def test_labelling_service_run_success(mock_db_port, mock_runner_port, mock_settings):
    atoms = Atoms('H', info={'db_id': 1})
    mock_db_port.get_atoms_by_state.return_value = [atoms]
    mock_runner_port.run.return_value = SAMPLE_SUCCESS_OUTPUT
    service = LabellingService(mock_db_port, mock_runner_port, mock_settings)
    service.run()
    mock_db_port.get_atoms_by_state.assert_called_once_with('unlabelled')
    mock_runner_port.run.assert_called_once()
    mock_db_port.update_with_dft_results.assert_called_once()
    call_args = mock_db_port.update_with_dft_results.call_args[0]
    assert call_args[0] == 1
    assert isinstance(call_args[1], DFTResult)
