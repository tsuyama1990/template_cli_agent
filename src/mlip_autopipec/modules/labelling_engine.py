import numpy as np
from ase.atoms import Atoms

from mlip_autopipec.domain_models import DFTResult


class LabellingEngine:
    """
    A placeholder for the Quantum Espresso labelling engine.
    In a real implementation, this class would call `pw.x` and parse its output.
    """

    def execute(self, atoms: Atoms) -> tuple[DFTResult, bool]:
        """
        Simulates a DFT calculation.
        Returns a DFTResult with dummy data and a success status.
        """
        print(f"Simulating DFT calculation for structure with {len(atoms)} atoms...")

        # Dummy data that mimics a successful calculation for a stable crystal
        dummy_energy = -150.0
        dummy_forces = (np.random.rand(*atoms.positions.shape) * 0.001).tolist()
        dummy_stress = (np.random.rand(6) * 0.1).tolist()

        result = DFTResult(
            energy=dummy_energy,
            forces=dummy_forces,
            stress=dummy_stress,
        )

        print("DFT calculation simulation successful.")
        return result, True
