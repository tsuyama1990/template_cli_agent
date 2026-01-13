from mlip_autopipec.common.pydantic_models import FullConfig
from mlip_autopipec.generators import AlloyGenerator
from mlip_autopipec.interfaces import ISampler, IStructureGenerator
from mlip_autopipec.samplers import RandomSampler


def create_generator(config: FullConfig) -> IStructureGenerator:
    """
    Factory function to create a structure generator based on the config.

    For Cycle 1, this only supports the AlloyGenerator.
    """
    # A more advanced implementation would have logic to select the
    # generator based on a 'type' field in the config.
    return AlloyGenerator(config.system)


def create_sampler(config: FullConfig) -> ISampler:
    """
    Factory function to create a sampler based on the config.

    For Cycle 1, this only supports the RandomSampler.
    """
    # A more advanced implementation would select the sampler based on
    # a 'method' field in a SamplingConfig Pydantic model.
    return RandomSampler()
