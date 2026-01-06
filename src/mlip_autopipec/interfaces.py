from abc import ABC, abstractmethod
from typing import Protocol

from ase import Atoms

from mlip_autopipec.common.pydantic_models import FullConfig


class IStructureGenerator(Protocol):
    """Interface for a structure generator."""

    def generate(self) -> list[Atoms]:
        """Generate a list of atomic structures."""
        ...


class ISampler(Protocol):
    """Interface for a sampler."""

    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """Sample a subset of structures from a trajectory."""
        ...


class IExplorer(Protocol):
    """Interface for an explorer engine."""

    def explore(self, atoms: Atoms) -> list[Atoms]:
        """Explore the potential energy surface starting from a given structure."""
        ...


class IDatabaseManager(Protocol):
    """Interface for a database manager."""

    def connect(self, path: str) -> None:
        """Connect to the database."""
        ...

    def write_structures(self, structures: list[Atoms]) -> None:
        """Write a list of structures to the database."""
        ...


class IOrchestrator(ABC):
    """Abstract base class for the main pipeline orchestrator."""

    def __init__(self, config: FullConfig) -> None:
        self.config = config

    @abstractmethod
    def run(self) -> None:
        """Run the full data generation pipeline."""
        raise NotImplementedError
