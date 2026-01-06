import numpy as np
from ase import Atoms

from mlip_autopipec.common.pydantic_models import MDConfig


class MDEngine:
    """
    A simplified Molecular Dynamics engine for Cycle 1.

    This class simulates a short MD run by applying a random perturbation
    to the atomic positions of the input structure.
    """

    def __init__(self, config: MDConfig) -> None:
        self.config = config

    def explore(self, atoms: Atoms) -> list[Atoms]:
        """
        "Explores" the structure by creating a few slightly modified copies.
        """
        trajectory = []
        for _ in range(5):  # Create a short, dummy trajectory
            frame = atoms.copy()  # type: ignore[no-untyped-call]
            positions = frame.get_positions()
            noise = np.random.uniform(-0.2, 0.2, size=positions.shape)
            frame.set_positions(positions + noise)
            frame.wrap()
            trajectory.append(frame)
        return trajectory
