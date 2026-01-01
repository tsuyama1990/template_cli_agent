import subprocess
from unittest.mock import patch

import pytest
from ase import Atoms

from mlip_autopipec.data.models import DFTCompute, DFTResults
from mlip_autopipec.modules.c_labeling_engine import LabelingEngine


@pytest.fixture
def dft_config():
    """Fixture for a sample DFTCompute config."""
    return DFTCompute(
        code="quantum_espresso",
        command="mpirun -np 4 pw.x",
        pseudopotentials={'Si': 'Si.UPF'},
        ecutwfc=60.0,
        ecutrho=240.0,
        kpoints_density=3.0,
    )


@patch('os.makedirs')
@patch('builtins.open')
@patch('subprocess.run')
@patch('mlip_autopipec.modules.c_labeling_engine.create_qe_input_from_atoms')
@patch('mlip_autopipec.modules.c_labeling_engine.parse_qe_output')
def test_labeling_engine_success(
    mock_parse, mock_create, mock_run, mock_open, mock_makedirs, dft_config
):
    """Test the LabelingEngine's successful execution path."""
    mock_create.return_value = "dummy_qe_input"
    mock_parse.return_value = DFTResults(
        energy=-1360.5, forces=[[1, 1, 1]], stress=[1, 2, 3, 4, 5, 6]
    )
    mock_run.return_value = subprocess.CompletedProcess(args=[], returncode=0)

    structures_to_label = [(1, Atoms('Si'))]
    engine = LabelingEngine(dft_config)
    results = engine.execute(structures_to_label)

    mock_create.assert_called_once()
    mock_run.assert_called_once()
    mock_parse.assert_called_once()

    assert len(results) == 1
    assert results[0][0] == 1
    assert isinstance(results[0][1], DFTResults)


@patch('os.makedirs')
@patch('builtins.open')
@patch('subprocess.run')
@patch('mlip_autopipec.modules.c_labeling_engine.create_qe_input_from_atoms')
def test_labeling_engine_dft_failure(
    mock_create, mock_run, mock_open, mock_makedirs, dft_config
):
    """Test that the engine raises an exception on DFT failure."""
    mock_run.side_effect = subprocess.CalledProcessError(
        returncode=1, cmd="pw.x", stderr=b"QE error"
    )

    structures_to_label = [(1, Atoms('Si'))]
    engine = LabelingEngine(dft_config)

    with pytest.raises(subprocess.CalledProcessError):
        engine.execute(structures_to_label)


@patch('os.makedirs')
@patch('builtins.open')
@patch('subprocess.run')
@patch('mlip_autopipec.modules.c_labeling_engine.create_qe_input_from_atoms')
@patch('mlip_autopipec.modules.c_labeling_engine.parse_qe_output')
def test_labeling_engine_parsing_failure(
    mock_parse, mock_create, mock_run, mock_open, mock_makedirs, dft_config
):
    """Test that the engine returns an empty list for parsing failures."""
    mock_parse.return_value = None
    mock_run.return_value = subprocess.CompletedProcess(args=[], returncode=0)

    structures_to_label = [(1, Atoms('Si'))]
    engine = LabelingEngine(dft_config)
    results = engine.execute(structures_to_label)

    mock_parse.assert_called_once()
    assert len(results) == 0
