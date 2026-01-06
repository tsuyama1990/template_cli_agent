# -*- coding: utf-8 -*-
"""Unit tests for the structure generators."""
from __future__ import annotations

from typing import Any

from mlip_autopipec.config.models import SystemConfig
from mlip_autopipec.generators.alloy import AlloyGenerator


def test_alloy_generator_produces_correct_structures() -> None:
    """Test that the AlloyGenerator creates structures with the correct properties."""
    # 1. Define a valid system configuration
    config_dict: dict[str, Any] = {
        "elements": ["Cu", "Au"],
        "composition": {"Cu": 0.5, "Au": 0.5},
        "lattice": "fcc",
        "num_structures": 5,
    }
    system_config = SystemConfig(**config_dict)

    # 2. Initialize and run the generator
    generator = AlloyGenerator(config=system_config)
    structures = generator.generate()

    # 3. Assert the output
    # Check that the correct number of structures was generated
    assert len(structures) == system_config.num_structures

    # Check the properties of each generated structure
    for atoms in structures:
        # Check composition
        symbols = atoms.get_chemical_symbols()  # type: ignore
        num_atoms = len(symbols)
        cu_count = symbols.count("Cu")
        au_count = symbols.count("Au")

        expected_cu = round(system_config.composition["Cu"] * num_atoms)
        expected_au = round(system_config.composition["Au"] * num_atoms)

        # Allow for small rounding differences in composition
        assert abs(cu_count - expected_cu) <= 1
        assert abs(au_count - expected_au) <= 1
        assert cu_count + au_count == num_atoms

        # Check for physical validity (no overlapping atoms)
        distances = atoms.get_all_distances(mic=True)  # type: ignore
        min_dist = min(d for d in distances.flatten() if d > 0)
        assert min_dist > 1.0, "Generated structure has overlapping atoms."
