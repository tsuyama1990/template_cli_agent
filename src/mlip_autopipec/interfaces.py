from abc import ABC, abstractmethod

from ase import Atoms


class BaseStructureGenerator(ABC):
    @abstractmethod
    def generate(self) -> list[Atoms]:
        ...

class BaseSampler(ABC):
    @abstractmethod
    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        ...
