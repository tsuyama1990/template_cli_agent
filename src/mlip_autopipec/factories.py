# src/mlip_autopipec/factories.py
"""
Factory functions for creating various components of the pipeline.
This pattern decouples the main orchestrator from the concrete implementations,
making it easy to add new generators or samplers in the future.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

from .generators.alloy import AlloyGenerator
from .interfaces import ISampler, IStructureGenerator
from .samplers.random_sampler import RandomSampler

if TYPE_CHECKING:
    from .common.pydantic_models import FullConfig


def create_generator(config: FullConfig) -> IStructureGenerator:
    """
    Factory function to create a structure generator instance.
    """
    return AlloyGenerator(config)


def create_sampler(config: FullConfig) -> ISampler:
    """
    Factory function to create a sampler instance.
    """
    method = config.sampling.method
    if method == "random":
        return RandomSampler(config)
    error_msg = f"Unsupported sampling method: {method}"
    raise ValueError(error_msg)
