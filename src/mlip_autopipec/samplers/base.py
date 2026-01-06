# src/mlip_autopipec/samplers/base.py
"""
Defines the abstract base class for all samplers.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from ase import Atoms

if TYPE_CHECKING:
    from mlip_autopipec.common.pydantic_models import FullConfig


class BaseSampler(ABC):
    """Abstract base class for sampling algorithms."""

    def __init__(self, config: FullConfig) -> None:
        self.config = config

    @abstractmethod
    def sample(self, trajectory: list[Atoms]) -> list[Atoms]:
        """Sample structures from a trajectory."""
        raise NotImplementedError
