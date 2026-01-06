# src/mlip_autopipec/explorers/md_engine.py
"""Molecular Dynamics and hybrid MD/MC simulation engine."""

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.common.pydantic_models import MDConfig


class MDEngine:
    """
    Manages the MD simulation part of the exploration stage.

    For Cycle 1, this will be a placeholder and will not run a real MD
    simulation. It will simply return the initial structures.
    """

    def __init__(self, config: MDConfig) -> None:
        self.config = config

    def run(self, initial_structures: list[Atoms]) -> list[Atoms]:
        """
        Runs the MD simulation.

        Args:
            initial_structures: A list of ASE Atoms objects to start from.

        Returns:
            A list of ASE Atoms objects representing the simulation trajectory.
            In Cycle 1, this is just the initial structures.
        """
        # In a real implementation, this would run a complex MD simulation.
        # For now, we just return a slightly modified copy of the input.
        trajectory = []
        for atoms in initial_structures:
            for _ in range(self.config.n_steps):
                # Simulate step by just rattling the atoms a bit more
                new_atoms = atoms.copy()  # type: ignore[no-untyped-call]
                new_atoms.rattle(stdev=0.05)  # type: ignore[no-untyped-call]
                # Attach a calculator to store the energy correctly
                calc = SinglePointCalculator(new_atoms, energy=-1.0)
                new_atoms.calc = calc
                trajectory.append(new_atoms)
        return trajectory
