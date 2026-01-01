import subprocess
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest
import yaml
from ase import Atoms
from click.testing import CliRunner

from mlip_autopipec.cli import app
from mlip_autopipec.database import AseDBWrapper
from mlip_autopipec.interfaces import IProcessRunner

DUMMY_QE_OUTPUT = """
     Program PWSCF v.6.5 starts on 12Jan2021 at 12:00:00
     lattice parameter (alat)  =  18.897260  a.u.
     celldm(1)   18.897260
     number of atoms/cell      =           3
     number of atomic types    =           2
     crystal axes: (angstrom)
               a(1) = (   10.000   0.000   0.000 )
               a(2) = (    0.000  10.000   0.000 )
               a(3) = (    0.000   0.000  10.000 )
     site n.     atom                  positions (alat units)
         1           O      tau(   1) = (   0.00000   0.00000   0.01170  )
         2           H      tau(   2) = (   0.00000   0.07570  -0.04690  )
         3           H      tau(   3) = (   0.00000  -0.07570  -0.04690  )
     Forces acting on atoms (cartesian axes, Ry/au):
          atom    1   type  O   force =     0.001   0.002   0.003
          atom    2   type  H   force =    -0.001  -0.001  -0.001
          atom    3   type  H   force =     0.000   0.000  -0.002
     total energy              =     -17.83123456 Ry
     total stress  (Ry/bohr**3)                (kbar)     P=       -0.01
      -0.00000011   0.00000000   0.00000000     -0.01      0.00      0.00
       0.00000000  -0.00000015   0.00000000      0.00     -0.02      0.00
       0.00000000   0.00000000  -0.00000020      0.00      0.00     -0.03
"""


@pytest.fixture
def mock_process_runner(mocker):
    """Mocks IProcessRunner and the ASE parser to simulate a QE execution."""

    # 1. Mock the subprocess call
    def mock_run(command, stdout_path):
        # The content doesn't matter, but the file must exist
        with open(stdout_path, "w") as f:
            f.write("dummy content")
        return subprocess.CompletedProcess(args=command, returncode=0)

    mock_runner = MagicMock(spec=IProcessRunner)
    mock_runner.run.side_effect = mock_run
    mocker.patch("mlip_autopipec.factories.SubprocessRunner", return_value=mock_runner)

    # 2. Mock the ASE parser
    mock_atoms = Atoms(
        "H2O",
        positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)],
        cell=np.eye(3) * 10,
        pbc=True,
    )
    # Energy in Ry, converted to eV for ASE
    energy = -17.83123456 * 13.6057
    forces = np.random.rand(3, 3)
    # ASE expects stress in a 6-element Voigt vector
    stress = np.random.rand(6)

    from ase.calculators.singlepoint import SinglePointCalculator

    mock_atoms.calc = SinglePointCalculator(
        atoms=mock_atoms, energy=energy, forces=forces, stress=stress
    )

    mocker.patch("mlip_autopipec.modules.labeling_engine.read", return_value=mock_atoms)

    return mock_runner


@pytest.fixture
def temp_alloy_config_file(tmp_path):
    """Creates a temporary user input config file for an alloy."""
    config_content = {
        "system": {"elements": ["Fe", "Pt"], "composition": "FePt"},
    }
    config_path = tmp_path / "input_alloy.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_content, f)
    return config_path


def test_full_workflow_from_composition(
    mock_process_runner, temp_alloy_config_file, mocker
):
    """
    Tests the full workflow starting from just a chemical composition.
    It verifies that the StructureGenerator runs, populates the database,
    and then the orchestrator proceeds to label the new structures.
    """
    # We also mock the training engine to avoid NotImplementedError
    mocker.patch(
        "mlip_autopipec.modules.training_engine.TrainingEngine.train", return_value=None
    )

    runner = CliRunner()
    with runner.isolated_filesystem() as temp_dir:
        config_path = Path(temp_dir) / "input.yaml"
        db_path = Path(temp_dir) / "test.db"

        config_content = temp_alloy_config_file.read_text()
        with open(config_path, "w") as f:
            f.write(config_content)

        # The database starts empty
        db_wrapper = AseDBWrapper(db_path=str(db_path))
        assert db_wrapper.is_empty()

        # Run the full workflow
        result = runner.invoke(
            app,
            ["run", "--config", str(config_path), "--db-path", str(db_path)],
            catch_exceptions=False,
        )
        assert result.exit_code == 0
        assert "Successfully generated and saved" in result.output

        # 1. Verify Database State after the full run
        # The database should be populated, and all generated structures should now be labeled.
        assert not db_wrapper.is_empty()
        labeled_data = db_wrapper.get_all_labeled_atoms()
        assert len(labeled_data) == 10
        unlabeled_ids = db_wrapper.get_unlabeled_ids()
        assert len(unlabeled_ids) == 0

        # 2. Verify Labeling Process
        # The labeling process should have been called for each new structure
        assert mock_process_runner.run.call_count == 10
        assert "Labeling complete for structure ID:" in result.output

        # The database should now have labeled structures
        labeled_data = db_wrapper.get_all_labeled_atoms()
        assert len(labeled_data) == 10

        # 3. Verify Training
        assert "Training process finished." in result.output
