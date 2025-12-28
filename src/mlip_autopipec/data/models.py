from pydantic import BaseModel


class DFTResult(BaseModel):
    """
    A Pydantic model to store the results of a DFT calculation.
    """

    total_energy_ev: float
    forces: list[list[float]]
    stress: list[list[float]]
    was_successful: bool
    error_message: str | None = None


class TrainingConfig(BaseModel):
    """
    A Pydantic model for MLIP training configuration.
    """

    model_type: str
    learning_rate: float
    num_epochs: int
    r_cut: float
    delta_learn: bool
    baseline_potential: str
