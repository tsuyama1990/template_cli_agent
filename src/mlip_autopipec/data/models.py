
from pydantic import BaseModel


class DFTResult(BaseModel):
    total_energy_ev: float | None = None
    forces: list[list[float]] | None = None
    stress: list[list[float]] | None = None
    was_successful: bool
    error_message: str | None = None

class TrainingConfig(BaseModel):
    model_type: str
    learning_rate: float
    num_epochs: int
    r_cut: float
    delta_learn: bool
    baseline_potential: str
