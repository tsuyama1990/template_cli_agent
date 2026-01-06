from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.interfaces import BaseSampler, BaseStructureGenerator
from mlip_autopipec.samplers.random_sampler import RandomSampler


def create_generator(config: FullConfig) -> BaseStructureGenerator:
    """
    Factory function to create a structure generator based on the config.

    For Cycle 1, this only supports the 'AlloyGenerator'. In the future, this
    function would inspect the config to decide which generator to instantiate.
    """
    # In a future cycle, you might have:
    # if config.system.type == 'alloy':
    #     return AlloyGenerator(config)
    # elif config.system.type == 'ionic':
    #     return IonicGenerator(config)
    # For now, we default to the only implemented generator.
    return AlloyGenerator(config)

def create_sampler(config: FullConfig) -> BaseSampler:
    """
    Factory function to create a sampler based on the config.
    """
    if config.sampling.method == "random":
        return RandomSampler(config)
    if config.sampling.method == "fps":
        # Placeholder for Farthest Point Sampling in Cycle 2
        raise NotImplementedError("FPS sampler is not yet implemented.")
    raise ValueError(f"Unknown sampling method: {config.sampling.method}")
