from typing import Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    model_validator,
)

ERROR_ELEMENTS_MISMATCH = "Elements and composition keys do not match."
ERROR_COMPOSITION_SUM = "Composition must sum to 1.0."


class SystemConfig(BaseModel):
    """
    Configuration for the physical system.

    Defines the material's composition, crystal structure, and the number of initial
    configurations to generate.
    """

    elements: list[str] = Field(
        ..., min_length=1, description="A list of chemical symbols (e.g., ['Fe', 'Pt'])."
    )
    composition: dict[str, float] = Field(
        ...,
        description="A mapping from element to its fractional composition (e.g., {'Fe': 0.5, 'Pt': 0.5}).",
    )
    lattice: Literal["fcc", "bcc", "hcp"] = Field(..., description="The crystal lattice type.")
    num_structures: int = Field(
        ..., gt=0, description="The number of initial seed structures to generate."
    )

    @model_validator(mode="after")
    def validate_system(self) -> "SystemConfig":
        """
        Validate the system configuration.

        Ensures that the elements listed match the keys in the composition dictionary
        and that the composition values sum to 1.0.
        """
        if set(self.elements) != set(self.composition.keys()):
            raise ValueError(ERROR_ELEMENTS_MISMATCH)
        if not abs(sum(self.composition.values()) - 1.0) < 1e-6:
            raise ValueError(ERROR_COMPOSITION_SUM)
        return self


class ExplorationConfig(BaseModel):
    """
    Configuration for the exploration stage.

    Contains settings for the simulation used to explore the potential energy surface.
    """

    temperature: float = Field(..., gt=0, description="The simulation temperature in Kelvin.")


class SamplingConfig(BaseModel):
    """
    Configuration for the sampling stage.

    Defines the method and parameters for selecting a subset of structures from the
    exploration stage.
    """

    method: Literal["random"] = Field(..., description="The sampling method to use.")
    fraction: float = Field(
        ...,
        gt=0,
        le=1,
        description="The fraction of structures to select from the exploration stage.",
    )


class FullConfig(BaseModel):
    """
    Top-level configuration model for the entire pipeline.

    This model aggregates all other configuration components and serves as the single
    source of truth for a pipeline run.
    """

    model_config = ConfigDict(extra="forbid")

    system: SystemConfig = Field(..., description="Defines the physical system to be generated.")
    exploration: ExplorationConfig = Field(
        ..., description="Defines the parameters for the exploration stage."
    )
    sampling: SamplingConfig = Field(
        ..., description="Defines the parameters for the sampling stage."
    )
    project_name: str = Field(
        ..., description="A unique name for the project or run, used for the output database file."
    )
