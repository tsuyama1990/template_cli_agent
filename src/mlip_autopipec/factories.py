from mlip_autopipec.common.pydantic_models import FullConfig, SamplingMethod
from mlip_autopipec.explorers.md_engine import MDEngine
from mlip_autopipec.generators.alloy import AlloyGenerator
from mlip_autopipec.interfaces import IExplorer, ISampler, IStructureGenerator
from mlip_autopipec.samplers.random_sampler import RandomSampler


def create_generator(config: FullConfig) -> IStructureGenerator:
    """Factory function to create a structure generator."""
    # In the future, this could have logic to choose different generators.
    return AlloyGenerator(
        system_config=config.system, generation_config=config.generation
    )


def create_explorer(config: FullConfig) -> IExplorer:
    """Factory function to create an explorer."""
    # In the future, this could have logic to choose different explorers.
    return MDEngine(config=config.exploration)


def create_sampler(config: FullConfig) -> ISampler:
    """Factory function to create a sampler."""
    if config.sampling.method == SamplingMethod.random:
        return RandomSampler(config=config.sampling)
    raise NotImplementedError(f"Sampler '{config.sampling.method}' not implemented.")
