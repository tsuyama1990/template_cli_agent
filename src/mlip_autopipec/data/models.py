"""Data models for the MLIP-AutoPipe project."""
from pydantic import BaseModel, Field
from typing import List, Optional

class DFTResult(BaseModel):
    """Represents the result of a DFT calculation."""
    total_energy_ev: float = Field(..., description="Total energy in eV")
    forces: List[List[float]] = Field(..., description="Forces on each atom in eV/Å")
    stress: List[List[float]] = Field(..., description="Stress tensor in eV/Å^3")
    was_successful: bool = Field(..., description="Whether the calculation was successful")
    error_message: Optional[str] = Field(None, description="Error message if the calculation failed")

class TrainingConfig(BaseModel):
    """Configuration for the MLIP training process."""
    model_type: str = Field(..., description="Type of the MLIP model (e.g., 'mace')")
    learning_rate: float = Field(..., description="Learning rate for the optimizer")
    num_epochs: int = Field(..., description="Number of training epochs")
    r_cut: float = Field(..., description="Cutoff radius for the potential")
    delta_learn: bool = Field(..., description="Whether to use Delta Learning")
    baseline_potential: str = Field(..., description="Name of the baseline potential for Delta Learning")
