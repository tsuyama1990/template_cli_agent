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

print("--- Running UAT-C01-02: Graceful DFT Failure Handling ---")

# GIVEN a physically problematic Atoms object
atoms = Atoms("Ar2", positions=[[0, 0, 0], [0, 0, 0.5]]) # Unphysically close

# AND a working QE installation is expected in the PATH
PSEUDO_DIR = os.environ.get("QE_PSEUDO_DIR", "./pseudos")
if not os.path.exists(PSEUDO_DIR):
    os.makedirs(PSEUDO_DIR)
    print(f"Created pseudo dir at {PSEUDO_DIR}")
    print("Please download SSSP pseudopotentials for Ar and place them there.")
    # A placeholder file
    with open(os.path.join(PSEUDO_DIR, "Ar.upf"), "w") as f:
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
print("1. Check that the script ran without crashing.")
print("2. Check that the log shows 'Labelling failed'.")
print("3. A 'model_cycle01.pt' file should NOT be created.")
print("4. An 'mlip.db' file should be created.")
print("   - Inspect it to confirm the structure state is 'labelling_failed'.")

assert not os.path.exists("model_cycle01.pt"), "UAT Failed: Model file was created."
assert os.path.exists("mlip.db"), "UAT Failed: Database file not created."

print("\n--- UAT-C01-02 Finished Successfully ---")
