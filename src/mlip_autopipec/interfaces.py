from abc import ABC, abstractmethod

from ase import Atoms


class IStructureGenerator(ABC):
    @abstractmethod
    def generate(self) -> list[Atoms]:
        ...


class IExplorer(ABC):
    @abstractmethod
    def run(self, initial_structures: list[Atoms]) -> list[Atoms]:
        ...


class ISampler(ABC):
    @abstractmethod
    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        ...
