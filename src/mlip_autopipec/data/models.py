from typing import Any, Dict, Literal

import numpy as np
from pydantic import BaseModel, Field, ConfigDict, field_validator


def ndarray_to_list(arr: np.ndarray) -> list:
    """Helper function to convert numpy arrays to lists for JSON serialization."""
    return arr.tolist()


class DFTResult(BaseModel):
    """
    A Pydantic model to store the results of a DFT calculation.
    Provides type safety and validation for DFT data.
    """

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        json_encoders={np.ndarray: ndarray_to_list},
    )

    energy: float = Field(..., description="The total potential energy calculated by DFT.")
    forces: np.ndarray = Field(..., description="Forces on each atom, shape (n_atoms, 3).")
    stress: np.ndarray = Field(..., description="Voigt-ordered stress tensor, shape (6,).")
    status: Literal["success", "failed"] = Field(..., description="The outcome of the calculation.")
    metadata: Dict = Field(default_factory=dict, description="Supplementary information.")

    @field_validator("forces", "stress", mode="before")
    @classmethod
    def validate_ndarray(cls, v: Any) -> np.ndarray:
        if isinstance(v, list):
            return np.array(v)
        return v
