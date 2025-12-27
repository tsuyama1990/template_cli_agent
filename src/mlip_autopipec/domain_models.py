from pydantic import BaseModel, ConfigDict


class DFTResult(BaseModel):
    """
    Represents the results of a DFT calculation.
    """

    model_config = ConfigDict(extra="forbid")

    energy: float
    forces: list[list[float]]
    stress: list[float] | list[list[float]] | None = None


class TrainingConfig(BaseModel):
    """
    Configuration for the MACE model training process.
    """

    model_config = ConfigDict(extra="forbid")

    model_path: str
    data_path: str
    epochs: int
    delta_learn: bool


class MainConfig(BaseModel):
    """
    Main configuration model for the application.
    """

    model_config = ConfigDict(extra="forbid")

    training: TrainingConfig
