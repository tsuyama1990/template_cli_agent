"""
Factories for creating concrete instances of components.

This module contains functions that act as factories to instantiate the
correct concrete implementation of a component based on the provided
configuration. This is a key part of the dependency injection and modular
design of the application.
"""

from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.interfaces import IGenerator, ISampler
from mlip_autopipec.samplers.random import RandomSampler


def create_generator(config: FullConfig) -> IGenerator:
    """
    Creates a structure generator based on the system configuration.

    Args:
        config: The full application configuration.

    Returns:
        An instance of a class that implements the IGenerator interface.

    Raises:
        ValueError: If an unsupported generator type is specified.
    """
    # For Cycle 1, we only have an alloy generator.
    # This factory can be extended to support other types in the future.
    # The type would be specified in the config, e.g., config.system.type
    generator_type = "alloy"  # Inferred for now

    if generator_type == "alloy":
        return AlloyGenerator(config.system)

    msg = f"Unsupported generator type: {generator_type}"
    raise ValueError(msg)


def create_sampler(config: FullConfig) -> ISampler:
    """
    Creates a sampler based on the sampling configuration.

    Args:
        config: The full application configuration.

    Returns:
        An instance of a class that implements the ISampler interface.

    Raises:
        ValueError: If an unsupported sampling method is specified.
    """
    method = config.sampling.method

    if method == "random":
        return RandomSampler(config.sampling)
    if method == "fps":
        msg = "Farthest Point Sampling (FPS) is a Cycle 2 feature."
        raise NotImplementedError(msg)

    msg = f"Unsupported sampling method: {method}"
    raise ValueError(msg)
