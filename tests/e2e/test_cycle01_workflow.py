import pytest
import os
import yaml
from ase import Atoms
from click.testing import CliRunner
from mlip_autopipec.main_cycle01 import main

# This test requires a Quantum Espresso installation and will be skipped if `pw.x` is not found.
# It also requires a simple pseudo potential file for Hydrogen.
QE_INSTALLED = os.system("which pw.x > /dev/null") == 0

@pytest.mark.skipif(not QE_INSTALLED, reason="Quantum Espresso (pw.x) not found in PATH.")
def test_e2e_cycle01_workflow():
    runner = CliRunner()

    with runner.isolated_filesystem():
        # 1. Create a dummy pseudo potential file
        os.makedirs("pseudos")
        with open("pseudos/H.upf", "w") as f:
            f.write("This is a dummy pseudo file.")

        # 2. Create a simple structure file (H2 molecule)
        atoms = Atoms('H2', positions=[[0, 0, 0], [0.74, 0, 0]], cell=[5, 5, 5])
        atoms.write("h2.xyz")

        # 3. Create a config file
        config = {
            "database_path": "test.db",
            "qe_command": "pw.x",
            "dft_parameters": {
                "pseudo_dir": "pseudos",
                "ecutwfc": 30.0,
                "kpoints": [1, 1, 1, 0, 0, 0],
                "pseudopotentials": {"H": "H.upf"}
            },
            "training": {
                "model_type": "mace",
                "learning_rate": 0.01,
                "num_epochs": 2, # Minimal epochs for a quick test
                "r_cut": 4.0,
                "delta_learn": True,
                "baseline_potential": "zbl"
            }
        }
        with open("config.yaml", "w") as f:
            yaml.dump(config, f)

        # 4. Run the CLI
        result = runner.invoke(main, ["--config", "config.yaml", "--structure", "h2.xyz"])

        # 5. Assertions
        assert result.exit_code == 0
        assert "Workflow complete." in result.output
        assert os.path.exists("models/trained_model.pt")

        # Check the database for success
        from mlip_autopipec.data.database import AseDB
        db = AseDB("test.db")
        db_entry = db.read(1)
        assert db_entry['was_successful'] is True
