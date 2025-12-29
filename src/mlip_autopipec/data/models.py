# ruff: noqa: D101, D102, D103, D104, D105, D107

from pydantic import BaseModel


class DFTResult(BaseModel):
    total_energy_ev: float
    forces: list[list[float]]
    stress: list[list[float]]
    was_successful: bool
    error_message: str | None = None


class TrainingConfig(BaseModel):
    model_type: str
    learning_rate: float
    num_epochs: int
    r_cut: float
    delta_learn: bool
    baseline_potential: str
