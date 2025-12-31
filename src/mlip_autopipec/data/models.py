"""Data models for MLIP-AutoPipe."""

import numpy as np
from pydantic import BaseModel, Field


class DFTResult(BaseModel):
    """
    Represents the result of a DFT calculation.
    """

    energy: float = Field(..., description="Total energy from the DFT calculation.")
    forces: np.ndarray = Field(..., description="Forces on each atom.")
    stress: np.ndarray = Field(..., description="Stress tensor of the system.")

    class Config:
        arbitrary_types_allowed = True
