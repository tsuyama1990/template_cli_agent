from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from ase import Atoms


def detect_vacuum(
    atoms: Atoms,
    slab_threshold_factor: float = 0.8,
) -> str:
    """
    Detects the system type (bulk, slab, or molecule) based on atomic coordinates.

    This robust statistical method analyzes the spread of atomic positions relative
    to the unit cell to distinguish between different system types.

    Args:
        atoms: The atomic structure to analyze.
        slab_threshold_factor: A factor used to determine if a system is a slab.
                               If the atomic layer thickness along one axis is
                               less than this factor times the cell dimension,
                               it's considered a slab.

    Returns:
        A string classification: "bulk", "slab", or "molecule".
    """
    # A non-periodic system is always a 'molecule' in a box.
    if not np.all(atoms.pbc):
        return "molecule"

    scaled_pos = atoms.get_scaled_positions()

    # Calculate the thickness of the atomic layer along each periodic direction.
    # This correctly handles coordinates that wrap around the periodic boundary.
    ranges = []
    for i in range(3):
        if atoms.pbc[i]:
            pos_1d = np.sort(scaled_pos[:, i])
            # The largest gap in coordinates indicates the vacuum thickness.
            gaps = np.diff(pos_1d)
            # Compare internal gaps with the gap across the boundary (1.0 - range).
            largest_gap = 1.0 - (pos_1d[-1] - pos_1d[0])
            if gaps.size > 0:
                largest_gap = max(largest_gap, np.max(gaps))
            # The atomic range is 1.0 minus the vacuum thickness.
            ranges.append(1.0 - largest_gap)
        else:
            ranges.append(np.ptp(scaled_pos[:, i]))

    # Classify based on how many axes have a large vacuum gap.
    is_slab_axis = [r < slab_threshold_factor for r in ranges]

    if sum(is_slab_axis) == 1:
        return "slab"  # One vacuum layer.
    if sum(is_slab_axis) > 1:
        return "molecule"  # e.g., a wire, or a cluster isolated on all sides.
    # sum(is_slab_axis) == 0
    return "bulk"  # Atoms are distributed throughout the cell.
