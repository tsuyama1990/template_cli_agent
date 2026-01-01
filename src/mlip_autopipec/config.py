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
)

# --- User Input Schemas (Minimal) ---


class UserInputSystem(BaseModel):
    """User-provided system information."""

    model_config = ConfigDict(extra="forbid")
    elements: list[str] = Field(..., min_length=1)
    composition: str


class UserInputSimulation(BaseModel):
    """User-provided simulation parameters (optional)."""

    model_config = ConfigDict(extra="forbid")
    temperature: list[int] | None = None


class UserInputConfig(BaseModel):
    """Top-level schema for the user's input.yaml."""

    model_config = ConfigDict(extra="forbid")
    system: UserInputSystem
    simulation: UserInputSimulation = UserInputSimulation()

    @classmethod
    def from_yaml(cls, path: str) -> "UserInputConfig":
        """Loads the configuration from a YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls(**data)


# --- Full Execution Schemas (Expanded) ---


class SystemConfig(BaseModel):
    """Expanded system configuration."""

    model_config = ConfigDict(extra="forbid")
    elements: list[str]
    composition: str
    structure_type: str
    melting_point_guess: float


class SimulationConfig(BaseModel):
    """Expanded simulation configuration."""

    model_config = ConfigDict(extra="forbid")
    temperature_steps: list[int]


class DFTComputeConfig(BaseModel):
    """Expanded DFT compute configuration."""

    model_config = ConfigDict(extra="forbid")
    ecutwfc: float
    ecutrho: float
    kpoints_density: float
    magnetism: str | None = None
    control: dict[str, Any]
    pseudopotentials: dict[str, str]


class ModelType(str, Enum):
    ACE = "ACE"


class MLIPTrainingConfig(BaseModel):
    """Pydantic model for MLIP training configuration."""

    model_config = ConfigDict(extra="forbid")
    model_type: ModelType
    r_cut: float
    loss_weights: dict[str, float]


class FullConfig(BaseModel):
    """Top-level schema for the fully expanded execution configuration."""

    model_config = ConfigDict(extra="forbid")
    system: SystemConfig
    simulation: SimulationConfig
    dft_compute: DFTComputeConfig
    mlip_training: MLIPTrainingConfig
    qe_command: str = "pw.x"
    db_path: str = "asedb.db"

    def to_yaml(self, path: str) -> None:
        """Saves the configuration to a YAML file."""
        with open(path, "w") as f:
            # Use mode='json' to ensure enums are serialized as strings
            yaml.dump(self.model_dump(mode="json"), f, sort_keys=False)


# --- Heuristic Engine ---

SSSP_PRECISION_CUTOFFS = {
    "H": (50, 400),
    "He": (60, 480),
    "Li": (40, 320),
    "Be": (50, 400),
    "B": (60, 480),
    "C": (60, 480),
    "N": (70, 560),
    "O": (70, 560),
    "F": (80, 640),
    "Ne": (90, 720),
    "Na": (40, 320),
    "Mg": (40, 320),
    "Al": (40, 320),
    "Si": (40, 320),
    "P": (50, 400),
    "S": (50, 400),
    "Cl": (60, 480),
    "Ar": (70, 560),
    "K": (40, 320),
    "Ca": (40, 320),
    "Sc": (60, 480),
    "Ti": (70, 560),
    "V": (70, 560),
    "Cr": (70, 560),
    "Mn": (70, 560),
    "Fe": (80, 640),
    "Co": (80, 640),
    "Ni": (80, 640),
    "Cu": (70, 560),
    "Zn": (70, 560),
    "Ga": (60, 480),
    "Ge": (60, 480),
    "As": (60, 480),
    "Se": (60, 480),
    "Br": (60, 480),
    "Kr": (70, 560),
    "Rb": (40, 320),
    "Sr": (40, 320),
    "Y": (60, 480),
    "Zr": (60, 480),
    "Nb": (70, 560),
    "Mo": (70, 560),
    "Tc": (70, 560),
    "Ru": (70, 560),
    "Rh": (70, 560),
    "Pd": (70, 560),
    "Ag": (60, 480),
    "Cd": (60, 480),
    "In": (60, 480),
    "Sn": (60, 480),
    "Sb": (60, 480),
    "Te": (60, 480),
    "I": (60, 480),
    "Xe": (70, 560),
    "Cs": (40, 320),
    "Ba": (40, 320),
    "La": (60, 480),
    "Ce": (60, 480),
    "Pr": (60, 480),
    "Nd": (60, 480),
    "Pm": (60, 480),
    "Sm": (60, 480),
    "Eu": (60, 480),
    "Gd": (70, 560),
    "Tb": (70, 560),
    "Dy": (70, 560),
    "Ho": (70, 560),
    "Er": (70, 560),
    "Tm": (70, 560),
    "Yb": (70, 560),
    "Lu": (70, 560),
    "Hf": (70, 560),
    "Ta": (70, 560),
    "W": (80, 640),
    "Re": (80, 640),
    "Os": (80, 640),
    "Ir": (80, 640),
    "Pt": (90, 720),
    "Au": (80, 640),
    "Hg": (70, 560),
    "Tl": (60, 480),
    "Pb": (60, 480),
    "Bi": (60, 480),
    "Po": (60, 480),
    "At": (60, 480),
    "Rn": (70, 560),
}

MAGNETIC_ELEMENTS = {"Fe", "Co", "Ni"}

ELEMENT_MELTING_POINTS = {
    "Fe": 1811,
    "Pt": 2041,
    "Si": 1687,
    "Al": 933,
}


class ConfigExpander:
    """Expands a minimal UserInputConfig into a comprehensive FullConfig."""

    def expand(self, user_input: UserInputConfig) -> FullConfig:
        """Applies heuristics to expand the user configuration."""

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
        if sim_input.temperature:
            steps = np.linspace(
                sim_input.temperature[0], sim_input.temperature[1], 3, dtype=int
            ).tolist()
        else:
            upper_bound = int(0.8 * melting_point)
            steps = np.linspace(300, upper_bound, 3, dtype=int).tolist()

        return SimulationConfig(temperature_steps=steps)

    def _expand_dft_config(self, system_input: UserInputSystem) -> DFTComputeConfig:
        max_ecutwfc = max(SSSP_PRECISION_CUTOFFS.get(el, (0, 0))[0] for el in system_input.elements)
        max_ecutrho = max(SSSP_PRECISION_CUTOFFS.get(el, (0, 0))[1] for el in system_input.elements)

        magnetism = (
            "ferromagnetic"
            if any(el in MAGNETIC_ELEMENTS for el in system_input.elements)
            else None
        )

        pseudos = {el: f"{el}.UPF" for el in system_input.elements}

        return DFTComputeConfig(
            ecutwfc=float(max_ecutwfc),
            ecutrho=float(max_ecutrho),
            kpoints_density=5.0,
            magnetism=magnetism,
            control={"calculation": "scf"},
            pseudopotentials=pseudos,
        )

    def _default_mlip_config(self) -> MLIPTrainingConfig:
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
    def serialize_numpy_array(self, v: np.ndarray) -> list:
        return v.tolist()
