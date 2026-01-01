"""
This module implements the structure generator for molecules.

It uses a Normal Mode Sampling (NMS) approach to generate a diverse set of
distorted structures, which is essential for exploring the potential energy
surface around an equilibrium geometry.
"""

import logging

import numpy as np
from ase import Atoms
from ase.build import molecule

from .base_generator import BaseStructureGenerator

logger = logging.getLogger(__name__)


class MoleculeGenerator(BaseStructureGenerator):
    """
    Generates initial structures for molecules using Normal Mode Sampling.
    """

    def generate(self) -> list[Atoms]:
        """
        Main method to generate a list of distorted molecular structures.

        This implementation is a simplified placeholder. A full implementation
        would require an external calculator to compute the Hessian matrix.
        Here, we simulate the process using random displacements as a proxy for
        true normal modes.

        Returns:
            A list of ase.Atoms objects representing distorted geometries.
        """
        base_atoms = self._get_equilibrium_geometry()

        # In a real scenario, we would calculate the Hessian and normal modes.
        # For this implementation, we simulate displacements.
        distorted_structures = self._apply_displacements(base_atoms)

        return [base_atoms] + distorted_structures

    def _get_equilibrium_geometry(self) -> Atoms:
        """
        Fetches or computes the equilibrium geometry of the molecule.

        For simplicity, this uses ase.build.molecule. A real implementation
        might perform a geometry optimization or query a database.

        Returns:
            An ase.Atoms object of the molecule at its equilibrium geometry.
        """
        try:
            atoms = molecule(self.system_info.composition)
            atoms.center(vacuum=3.0)
            logger.info(f"Created equilibrium geometry for {self.system_info.composition}")
            return atoms
        except Exception as e:
            raise ValueError(
                f"Could not create molecule '{self.system_info.composition}'. "
                "Ensure it is a valid chemical formula known to ASE."
            ) from e

    def _apply_displacements(self, base_atoms: Atoms) -> list[Atoms]:
        """
        Applies random displacements to the atoms to simulate normal modes.
        """
        distorted_structures = []
        num_structures = self.config.simulation.initial_structures_to_generate - 1

        for _ in range(num_structures):
            distorted_atoms = base_atoms.copy()
            # Apply random displacements with a magnitude of 0.1 Angstrom
            displacements = np.random.rand(*distorted_atoms.positions.shape) * 0.2 - 0.1
            distorted_atoms.positions += displacements
            distorted_structures.append(distorted_atoms)

        logger.info(f"Generated {num_structures} distorted structures via random displacement.")
        return distorted_structures
