"""
This module defines the Pydantic-based configuration system for the application.

It implements a two-tier configuration strategy:
1.  **UserInputConfig**: A minimal set of parameters provided by the user in a
    simple `input.yaml` file. This is designed to be easy for non-experts.
2.  **FullConfig**: A comprehensive, fully-populated configuration object that
    contains every parameter needed for the entire workflow to run.

The `ConfigExpander` class acts as a heuristic engine to transform the
`UserInputConfig` into the `FullConfig`, applying physics-based heuristics
and sensible defaults. This module is central to the project's goal of
"removing the human expert from the loop."
"""

from enum import Enum
from typing import Any

import numpy as np
import yaml
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    field_serializer,
    field_validator,
    model_validator,
)

from .constants import (
    ELEMENT_MELTING_POINTS,
    MAGNETIC_ELEMENTS,
    SSSP_PRECISION_CUTOFFS,
)

# --- User Input Schemas (Minimal) ---


class UserInputSystem(BaseModel):
    """
    Defines the chemical system as provided by the user.

    Attributes:
        elements: A list of chemical symbols (e.g., ['Fe', 'Pt']).
        composition: A string representing the stoichiometry (e.g., 'FePt').
    """

    model_config = ConfigDict(extra="forbid")
    elements: list[str] = Field(..., min_length=1)
    composition: str


class UserInputSimulation(BaseModel):
    """
    Optional simulation parameters provided by the user.

    Attributes:
        temperature: An optional list of two integers defining a temperature range.
    """

    model_config = ConfigDict(extra="forbid")
    temperature: list[int] | None = None


class UserInputConfig(BaseModel):
    """
    Top-level schema for parsing the user's `input.yaml` file.

    Attributes:
        system: The chemical system definition.
        simulation: Optional simulation parameters.
    """

    model_config = ConfigDict(extra="forbid")
    system: UserInputSystem
    simulation: UserInputSimulation = UserInputSimulation()

    @classmethod
    def from_yaml(cls, path: str) -> "UserInputConfig":
        """
        Loads and validates the user configuration from a YAML file.

        Args:
            path: The path to the `input.yaml` file.

        Returns:
            A validated UserInputConfig instance.
        """
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls(**data)


# --- Full Execution Schemas (Expanded) ---


class SystemConfig(BaseModel):
    """
    Expanded and detailed configuration for the chemical system.

    Attributes:
        elements: A list of chemical symbols.
        composition: The stoichiometry string.
        structure_type: The type of material (e.g., 'alloy', 'covalent').
        melting_point_guess: A heuristic-based guess for the melting point in Kelvin.
    """

    model_config = ConfigDict(extra="forbid")
    elements: list[str]
    composition: str
    structure_type: str
    melting_point_guess: float


class SimulationConfig(BaseModel):
    """
    Expanded configuration for simulation parameters.

    Attributes:
        temperature_steps: A list of temperatures (in Kelvin) for simulations.
        initial_structures_to_generate: The number of structures to generate.
    """

    model_config = ConfigDict(extra="forbid")
    temperature_steps: list[int]
    initial_structures_to_generate: int = 10


class DFTComputeConfig(BaseModel):
    """
    Detailed configuration for the DFT (Quantum Espresso) calculations.

    Attributes:
        ecutwfc: Plane-wave cutoff energy for wavefunctions (in Ry).
        ecutrho: Plane-wave cutoff energy for charge density (in Ry).
        kpoints_density: Density of k-points for Brillouin zone sampling.
        magnetism: Type of magnetism to simulate (e.g., 'ferromagnetic').
        control: A dictionary of Quantum Espresso control parameters.
        pseudopotentials: A mapping of element symbols to pseudopotential filenames.
    """

    model_config = ConfigDict(extra="forbid")
    ecutwfc: float
    ecutrho: float
    kpoints_density: float
    magnetism: str | None = None
    control: dict[str, Any]
    pseudopotentials: dict[str, str]


class ModelType(str, Enum):
    """Enumeration for the type of Machine Learned Interatomic Potential."""

    ACE = "ACE"


class MLIPTrainingConfig(BaseModel):
    """
    Configuration for training the MLIP model.

    Attributes:
        model_type: The type of MLIP model to train.
        r_cut: The cutoff radius for the potential.
        loss_weights: Weights for the energy and forces components of the loss function.
    """

    model_config = ConfigDict(extra="forbid")
    model_type: ModelType
    r_cut: float
    loss_weights: dict[str, float]


class FullConfig(BaseModel):
    """
    Top-level schema for the fully expanded and validated execution configuration.

    This object contains all parameters necessary to run the entire workflow and
    serves as the single source of truth for all components.

    Attributes:
        system: Detailed system configuration.
        simulation: Detailed simulation configuration.
        dft_compute: Detailed DFT calculation configuration.
        mlip_training: Detailed MLIP training configuration.
        qe_command: The command to execute Quantum Espresso.
        db_path: The path to the ASE database file.
    """

    model_config = ConfigDict(extra="forbid")
    system: SystemConfig
    simulation: SimulationConfig
    dft_compute: DFTComputeConfig
    mlip_training: MLIPTrainingConfig
    qe_command: str = "pw.x"
    db_path: str = "asedb.db"

    @model_validator(mode="after")
    def validate_pseudos_exist_for_all_elements(self) -> "FullConfig":
        """
        Ensures that a pseudopotential is defined for every element in the system.
        This acts as a cross-field validation check for configuration integrity.
        """
        elements = set(self.system.elements)
        pseudo_elements = set(self.dft_compute.pseudopotentials.keys())
        if not elements.issubset(pseudo_elements):
            missing = elements - pseudo_elements
            raise ValueError(f"Missing pseudopotentials for elements: {sorted(list(missing))}")
        return self

    def to_yaml(self, path: str) -> None:
        """
        Saves the configuration to a YAML file for inspection and reproducibility.

        Args:
            path: The path to save the YAML file to.
        """
        with open(path, "w") as f:
            # Use mode='json' to ensure enums are serialized as strings
            yaml.dump(self.model_dump(mode="json"), f, sort_keys=False)


# --- Heuristic Engine ---


class ConfigExpander:
    """
    Expands a minimal UserInputConfig into a comprehensive FullConfig.

    This class contains the core heuristic logic of the application. It takes
    the user's simple input and applies a series of rules and defaults to
    generate a complete configuration suitable for execution.
    """

    def expand(self, user_input: UserInputConfig) -> FullConfig:
        """
        Applies heuristics to expand the user configuration.

        Args:
            user_input: The validated user input configuration.

        Returns:
            A fully expanded and validated configuration object.
        """

        system = self._expand_system_config(user_input.system)
        simulation = self._expand_simulation_config(
            user_input.simulation, system.melting_point_guess
        )
        dft_compute = self._expand_dft_config(user_input.system)
        mlip_training = self._default_mlip_config()

        return FullConfig(
            system=system,
            simulation=simulation,
            dft_compute=dft_compute,
            mlip_training=mlip_training,
        )

    def _expand_system_config(self, system_input: UserInputSystem) -> SystemConfig:
        """Generates a detailed system configuration from user input."""
        structure_type = "covalent" if len(system_input.elements) == 1 else "alloy"

        melting_point = np.mean(
            [ELEMENT_MELTING_POINTS.get(el, 1500) for el in system_input.elements]
        )

        return SystemConfig(
            elements=system_input.elements,
            composition=system_input.composition,
            structure_type=structure_type,
            melting_point_guess=melting_point,
        )

    def _expand_simulation_config(
        self, sim_input: UserInputSimulation, melting_point: float
    ) -> SimulationConfig:
        """Generates a detailed simulation configuration from user input."""
        if sim_input.temperature:
            # If the user provides a temperature range, interpolate it into steps.
            steps = np.linspace(
                sim_input.temperature[0], sim_input.temperature[1], 3, dtype=int
            ).tolist()
        else:
            # Otherwise, create a default range from 300K up to 80% of the
            # estimated melting point, which is a common heuristic for exploring
            # the solid phase of a material.
            upper_bound = int(0.8 * melting_point)
            steps = np.linspace(300, upper_bound, 3, dtype=int).tolist()

        return SimulationConfig(temperature_steps=steps)

    def _expand_dft_config(self, system_input: UserInputSystem) -> DFTComputeConfig:
        """Generates a detailed DFT configuration using heuristics."""
        try:
            # The SSSP standard provides recommended cutoffs for each element.
            # To ensure accuracy, we must use the maximum required cutoff
            # among all elements present in the system.
            max_ecutwfc = max(SSSP_PRECISION_CUTOFFS[el][0] for el in system_input.elements)
            max_ecutrho = max(SSSP_PRECISION_CUTOFFS[el][1] for el in system_input.elements)
        except KeyError as e:
            raise ValueError(
                f"Unknown element '{e.args[0]}' provided. Cannot determine SSSP cutoff."
            ) from e

        # Enable spin-polarized calculations automatically if known magnetic
        # elements are present.
        magnetism = (
            "ferromagnetic"
            if any(el in MAGNETIC_ELEMENTS for el in system_input.elements)
            else None
        )

        pseudos = {el: f"{el}.UPF" for el in system_input.elements}

        return DFTComputeConfig(
            ecutwfc=float(max_ecutwfc),
            ecutrho=float(max_ecutrho),
            # A k-point density of 5.0 is a reasonable default for many systems,
            # balancing accuracy and computational cost.
            kpoints_density=5.0,
            magnetism=magnetism,
            control={"calculation": "scf"},
            pseudopotentials=pseudos,
        )

    def _default_mlip_config(self) -> MLIPTrainingConfig:
        """Returns a default MLIP training configuration."""
        return MLIPTrainingConfig(
            model_type=ModelType.ACE,
            r_cut=5.0,
            loss_weights={"energy": 1.0, "forces": 100.0},
        )


# --- Legacy Models ---


class DFTResult(BaseModel):
    """Pydantic model for storing DFT calculation results."""

    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)
    energy: float = Field(..., description="Total energy from the calculation.")
    forces: np.ndarray = Field(..., description="Forces on each atom.")
    stress: np.ndarray = Field(..., description="Stress tensor.")

    @field_validator("forces", "stress", mode="before")
    @classmethod
    def validate_numpy_array(cls, v: Any) -> np.ndarray:
        if isinstance(v, list):
            return np.array(v)
        if not isinstance(v, np.ndarray):
            raise ValueError("Value must be a numpy array or a list.")
        return v

    @field_serializer("forces", "stress")
    def serialize_numpy_array(self, v: np.ndarray) -> list[Any]:
        return v.tolist()
