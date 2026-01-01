import json
import re
from pathlib import Path
from typing import Dict, Any

import yaml
from mlip_autopipec.data.models import (
    FullConfig,
    MinimalSystem,
    System,
    StructureGeneration,
    DFTCompute,
    MLIPTraining,
)


class ConfigExpander:
    """
    Expands a minimal user configuration into a full, explicit configuration
    by applying a set of heuristics and default values.
    """

    def __init__(self, defaults_path: Path):
        """
        Initializes the ConfigExpander with the path to the defaults file.

        Args:
            defaults_path: Path to the JSON file containing default values,
                           e.g., SSSP pseudopotential recommendations.
        """
        if not defaults_path.exists():
            raise FileNotFoundError(f"Defaults file not found at: {defaults_path}")
        with open(defaults_path) as f:
            self._defaults = json.load(f)

    def expand_config(self, input_path: Path) -> FullConfig:
        """
        Loads a minimal input YAML, applies heuristics, and returns a
        validated FullConfig object.

        Args:
            input_path: Path to the user's minimal `input.yaml` file.

        Returns:
            A fully populated and validated FullConfig object.
        """
        # 1. Load and validate the minimal input
        with open(input_path) as f:
            data = yaml.safe_load(f)
        minimal_system = MinimalSystem(**data["system"])

        # 2. Apply Heuristics
        # System section
        elements = sorted(minimal_system.elements)
        composition_dict = self._parse_composition(minimal_system.composition)
        system = System(
            elements=elements,
            composition=composition_dict,
            structure_type="alloy",  # Hardcoded for now as per spec
        )

        # Generation section (using default values)
        generation = StructureGeneration(
            generation_strategy="sqs",
            supercell_size=16,
            strains=[-0.02, -0.01, 0.0, 0.01, 0.02],
        )

        # DFT Compute section
        dft_compute = self._infer_dft_params(elements)

        # MLIP Training section (using default values)
        mlip_training = MLIPTraining(
            model_type="mace",
            r_cut=5.0,
            delta_learning=False,
            loss_weights={"energy": 1.0, "forces": 100.0},
        )

        # 3. Assemble and validate the FullConfig
        full_config = FullConfig(
            system=system,
            generation=generation,
            dft_compute=dft_compute,
            mlip_training=mlip_training,
        )

        return full_config

    def _parse_composition(self, composition: Any) -> Dict[str, int]:
        """Parses a composition string like "FePt" into a dictionary."""
        if isinstance(composition, dict):
            return composition
        # Use regex to find elements and their optional counts
        # e.g., "Fe2Pt" -> [('Fe', '2'), ('Pt', '')]
        parts = re.findall(r"([A-Z][a-z]*)(\d*)", composition)
        comp_dict = {}
        for element, count in parts:
            comp_dict[element] = int(count) if count else 1
        return comp_dict

    def _infer_dft_params(self, elements: list[str]) -> DFTCompute:
        """Infers DFT parameters based on the elements in the system."""
        missing_elements = [el for el in elements if el not in self._defaults]
        if missing_elements:
            raise ValueError(f"Missing defaults for elements: {missing_elements}")

        max_ecutwfc = max(self._defaults[el]["ecutwfc"] for el in elements)
        max_ecutrho = max(self._defaults[el]["ecutrho"] for el in elements)
        pseudos = {el: self._defaults[el]["filename"] for el in elements}

        return DFTCompute(
            code="quantum_espresso",
            command="pw.x",
            pseudopotentials=pseudos,
            ecutwfc=max_ecutwfc,
            ecutrho=max_ecutrho,
            kpoints_density=3.0,  # A sensible default
        )
