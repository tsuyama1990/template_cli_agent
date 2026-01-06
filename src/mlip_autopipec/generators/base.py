# -*- coding: utf-8 -*-
"""
Abstract BaseStructureGenerator class, defining the interface for all structure generators.

This module defines the abstract base class that all structure generators must
inherit from, ensuring a consistent interface across the application.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms


class BaseStructureGenerator(ABC):
    """
    Abstract base class for all structure generators.

    This class defines the contract for all structure generation components.
    Each concrete generator must implement the `generate` method.
    """

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """
        Generate a list of atomic structures.

        This method should be implemented by all concrete generator classes.

        Returns:
            A list of ase.Atoms objects representing the generated structures.
        """
        raise NotImplementedError
