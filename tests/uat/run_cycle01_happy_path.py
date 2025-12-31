import os
import shutil

from ase import Atoms
from mlip_autopipec.orchestrator import Orchestrator

# Clean up previous runs
if os.path.exists("mlip.db"):
    os.remove("mlip.db")
if os.path.exists("model_cycle01.pt"):
    os.remove("model_cycle01.pt")
if os.path.exists("./out"):
    shutil.rmtree("./out")


print("--- Running UAT-C01-01: Successful End-to-End Workflow ---")

# GIVEN a simple, valid Atoms object
atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])

# AND a working QE installation is expected in the PATH
# User needs to provide the path to their pseudopotentials
PSEUDO_DIR = os.environ.get("QE_PSEUDO_DIR", "./pseudos")
if not os.path.exists(PSEUDO_DIR):
    os.makedirs(PSEUDO_DIR)
    print(f"Created pseudo dir at {PSEUDO_DIR}")
    print("Please download SSSP pseudopotentials for H and place them there.")
    # A placeholder file
    with open(os.path.join(PSEUDO_DIR, "H.upf"), "w") as f:
        f.write("This is a placeholder. Replace with a real pseudopotential file.")


orchestrator = Orchestrator(
    db_path="mlip.db",
    qe_command="pw.x", # Assumes pw.x is in the user's PATH
    pseudo_dir=PSEUDO_DIR,
)

# WHEN the workflow is executed
orchestrator.run_cycle01_workflow(atoms)

# THEN the user should verify the outputs
print("\n--- UAT Verification ---")
print("1. Check that the script ran without errors.")
print("2. A 'model_cycle01.pt' file should be created in this directory.")
print("3. An 'mlip.db' file should be created and contain the results.")
print("   - You can inspect it with a SQLite browser.")

assert os.path.exists("model_cycle01.pt"), "UAT Failed: Model file not created."
assert os.path.exists("mlip.db"), "UAT Failed: Database file not created."

print("\n--- UAT-C01-01 Finished Successfully ---")
