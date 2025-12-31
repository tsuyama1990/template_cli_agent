from typing import Any

import torch
from ase import Atoms
from ase.calculators.lj import LennardJones
from ase.calculators.singlepoint import SinglePointCalculator

from mlip_autopipec.config.models import TrainingParams


class TrainingEngine:
    """Encapsulates all model training logic."""

    def __init__(self, train_config: TrainingParams):
        """
        Initializes the TrainingEngine.

        Args:
            train_config: The training-specific configuration.
        """
        self.train_config = train_config

    def run(self, dataset: list[Atoms]):
        """
        Runs the training process on a list of labelled Atoms objects.

        Args:
            dataset: The labelled dataset for training.
        """
        formatted_data = self._prepare_data(dataset)
        trained_model = self._train_model(formatted_data)
        self._save_model(trained_model)

    def _prepare_data(self, dataset: list[Atoms]) -> Any:
        """
        Prepares the dataset for the training library.
        If delta_learning is enabled, it calculates the difference between
        the DFT labels and a baseline potential.
        """
        if not self.train_config.delta_learning:
            # In a real-world scenario, we would format this for MACE
            return dataset

        baseline_potential = LennardJones()
        delta_dataset = []

        for atoms in dataset:
            # Ensure the atoms object has a calculator with results
            if not atoms.calc:
                continue

            # Get the high-accuracy DFT results
            dft_energy = atoms.get_potential_energy()
            dft_forces = atoms.get_forces()

            # Calculate the low-accuracy baseline results
            atoms_copy = atoms.copy()
            atoms_copy.calc = baseline_potential
            baseline_energy = atoms_copy.get_potential_energy()
            baseline_forces = atoms_copy.get_forces()

            # Calculate the delta
            delta_energy = dft_energy - baseline_energy
            delta_forces = dft_forces - baseline_forces

            # Create a new Atoms object with the delta results
            delta_atoms = atoms.copy()
            delta_atoms.calc = SinglePointCalculator(
                delta_atoms, energy=delta_energy, forces=delta_forces
            )
            delta_dataset.append(delta_atoms)

        # In a real implementation, this list would be converted
        # to the format required by the MACE library. For the unit test,
        # returning the list of Atoms with delta calculators is sufficient.
        return delta_dataset

    def _train_model(self, formatted_data: Any) -> Any:
        """
        Wraps the underlying MACE/ACE library's fitting function.
        """
        # Placeholder: a mock model object
        print("Training model...")
        mock_model = torch.nn.Linear(10, 1)
        print("Training complete.")
        return mock_model

    def _save_model(self, trained_model: Any):
        """
        Serializes and saves the final trained potential to a file.
        """
        # Placeholder implementation
        output_path = "model.pt"
        print(f"Saving model to {output_path}...")
        torch.save(trained_model.state_dict(), output_path)
        print("Model saved.")
