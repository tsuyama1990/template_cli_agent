import pytest


@pytest.mark.skip(reason="Implementation not yet available")
def test_random_sampler_selects_correct_number():
    """Tests that the RandomSampler returns the requested number of samples."""
    # from mlip_autopipec.samplers.random_sampler import RandomSampler
    #
    # config = MagicMock()
    # config.sampling.num_samples = 5
    # sampler = RandomSampler(config)
    #
    # # Create a dummy trajectory with more structures than num_samples
    # dummy_trajectory = [Atoms('H', positions=[(0, 0, i)]) for i in range(10)]
    #
    # sampled_structures = sampler.sample(dummy_trajectory)
    #
    # assert len(sampled_structures) == 5
    # assert all(isinstance(s, Atoms) for s in sampled_structures)

@pytest.mark.skip(reason="Implementation not yet available")
def test_random_sampler_handles_small_trajectory():
    """Tests that the sampler returns all structures if the trajectory is smaller than num_samples."""
    # from mlip_autopipec.samplers.random_sampler import RandomSampler
    #
    # config = MagicMock()
    # config.sampling.num_samples = 10
    # sampler = RandomSampler(config)
    #
    # # Create a dummy trajectory with fewer structures than num_samples
    # dummy_trajectory = [Atoms('H', positions=[(0, 0, i)]) for i in range(5)]
    #
    # sampled_structures = sampler.sample(dummy_trajectory)
    #
    # assert len(sampled_structures) == 5
