import subprocess
from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pytest
import yaml
from ase import Atoms
from click.testing import CliRunner

from mlip_autopipec.cli import app
from mlip_autopipec.config import FullConfig
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
def temp_config_file(tmp_path):
    """Creates a temporary user input config file."""
    config_content = {
        "system": {"elements": ["H", "O"], "composition": "H2O"},
        "simulation": {"temperature": [300, 500]},
    }
    config_path = tmp_path / "input.yaml"
    with open(config_path, "w") as f:
        yaml.dump(config_content, f)
    return config_path


def test_label_and_train_workflow(mock_process_runner, temp_config_file):
    """An integration test for the full label-and-train workflow."""
    runner = CliRunner()
    with runner.isolated_filesystem() as temp_dir:
        config_path = Path(temp_dir) / "input.yaml"
        db_path = Path(temp_dir) / "test.db"

        config_content = temp_config_file.read_text()
        with open(config_path, "w") as f:
            f.write(config_content)

        db_wrapper = AseDBWrapper(db_path=str(db_path))
        atoms = Atoms("H2O", positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0)])
        structure_id = db_wrapper.add_atoms(atoms)
        assert structure_id == 1

        result = runner.invoke(
            app,
            [
                "label",
                "--config",
                str(config_path),
                "--id",
                str(structure_id),
                "--db-path",
                str(db_path),
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0
        assert f"Labeling complete for structure ID: {structure_id}" in result.output

        dump_path = Path(temp_dir) / "exec_config_dump.yaml"
        assert dump_path.exists()
        with open(dump_path) as f:
            dump_data = yaml.safe_load(f)
            parsed_config = FullConfig.model_validate(dump_data)
            assert parsed_config.db_path == str(db_path)

        mock_process_runner.run.assert_called_once()

        labeled_data = db_wrapper.get_all_labeled_atoms()
        assert len(labeled_data) == 1
        _, dft_result = labeled_data[0]

        # Verify the parsed data stored in the database
        expected_energy_ev = -17.83123456 * 13.6057
        assert dft_result.energy == pytest.approx(expected_energy_ev)
        assert dft_result.forces.shape == (3, 3)
        assert dft_result.stress.shape == (3, 3)

        result = runner.invoke(
            app,
            ["train", "--config", str(config_path), "--db-path", str(db_path)],
            catch_exceptions=False,
        )
        assert result.exit_code == 0
        assert "Training not implemented" in result.output
