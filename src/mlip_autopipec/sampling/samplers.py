# -*- coding: utf-8 -*-
"""
Implementation of sampling algorithms.

This module provides classes for various structure sampling strategies.
For Cycle 1, only a simple random sampler is implemented.
"""
from __future__ import annotations

import random
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ase import Atoms

    from mlip_autopipec.config.models import SamplingConfig


class RandomSampler:
    """
    A simple sampler that selects a random subset of structures.

    This sampler implements a basic random sampling strategy where a specified
    fraction of the input structures are chosen without replacement.
    """

    def __init__(self, config: SamplingConfig) -> None:
        """
        Initialize the RandomSampler.

        Args:
            config: The sampling configuration Pydantic model.
        """
        self.config = config

    def sample(self, structures: list[Atoms]) -> list[Atoms]:
        """
        Select a random subset of the provided structures.

        Args:
            structures: A list of ase.Atoms objects to be sampled from.

        Returns:
            A new list containing a randomly selected subset of the input structures.
        """
        num_to_select = int(len(structures) * self.config.fraction)
        return random.sample(structures, num_to_select)
