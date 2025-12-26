from typing import List
import torch
import numpy as np
from ase.atoms import Atoms
from mace.data.atomic_data import AtomicData
from mace.tools import AtomicNumberTable
from mace.data.neighborhood import get_neighborhood

def to_one_hot(tensor: torch.Tensor, num_classes: int) -> torch.Tensor:
    """Converts a tensor of indices to a one-hot tensor."""
    return torch.nn.functional.one_hot(tensor, num_classes=num_classes).to(torch.float32)

def to_atomic_data(
    atoms: Atoms,
    z_table: AtomicNumberTable,
    cutoff: float,
) -> AtomicData:
    """Converts an ASE Atoms object to a MACE AtomicData object."""

    # Get neighborhood info
    positions = atoms.get_positions()
    edge_index, shifts, unit_shifts, cell = get_neighborhood(
        positions=positions, # Pass numpy array directly
        cutoff=cutoff,
        pbc=atoms.pbc,
        cell=atoms.cell,
    )

    # Get node attributes (one-hot encoding of atomic numbers)
    atomic_numbers = atoms.get_atomic_numbers()
    indices = torch.tensor([z_table.z_to_index(z) for z in atomic_numbers], dtype=torch.long)
    node_attrs = to_one_hot(indices, num_classes=len(z_table))

    # Safely get energy and forces, providing defaults if not present
    energy = atoms.info.get('energy', 0.0)
    forces = atoms.arrays.get('forces', np.zeros((len(atoms), 3)))

    return AtomicData(
        edge_index=torch.tensor(edge_index, dtype=torch.long),
        node_attrs=node_attrs,
        positions=torch.tensor(positions, dtype=torch.float32),
        shifts=torch.tensor(shifts, dtype=torch.float32),
        unit_shifts=torch.tensor(unit_shifts, dtype=torch.float32),
        cell=torch.tensor(cell, dtype=torch.float32),
        weight=torch.tensor(1.0, dtype=torch.float32),
        forces=torch.tensor(forces, dtype=torch.float32),
        energy=torch.tensor(energy, dtype=torch.float32),
        stress=torch.zeros((1, 3, 3), dtype=torch.float32),
        virials=torch.zeros((1, 3, 3), dtype=torch.float32),
        dipole=torch.zeros((1, 3), dtype=torch.float32),
        charges=torch.zeros(len(atoms), dtype=torch.float32),
        head=None,
        energy_weight=None,
        forces_weight=None,
        stress_weight=None,
        virials_weight=None,
        dipole_weight=None,
        charges_weight=None,
        polarizability_weight=None,
        polarizability=None,
        elec_temp=None,
    )
