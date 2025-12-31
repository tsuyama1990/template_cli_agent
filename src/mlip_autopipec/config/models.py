
from pydantic import BaseModel, Field


class DFTParams(BaseModel):
    """Parameters for the DFT calculation using Quantum Espresso."""

    command: str = Field(
        ..., description="The command to execute Quantum Espresso, e.g., 'pw.x'."
    )
    pseudopotentials: dict[str, str] = Field(
        ..., description="Mapping of element symbol to pseudopotential filename."
    )
    ecutwfc: float = Field(
        ..., description="Wavefunction kinetic energy cutoff in Ry."
    )
    kpoints_density: float = Field(
        default=5000,
        description="K-point density for automatic k-point mesh generation.",
    )


class TrainingParams(BaseModel):
    """Parameters for the MLIP training."""

    model_type: str = Field(default="ace", description="The type of MLIP model to train.")
    r_cut: float = Field(..., description="Cutoff radius for the potential.")
    delta_learning: bool = Field(
        default=True, description="Whether to use delta learning."
    )


class FullConfig(BaseModel):
    """The full, validated configuration for the pipeline."""

    dft_compute: DFTParams
    training: TrainingParams
    ase_db_path: str = Field(..., description="Path to the ASE database file.")
