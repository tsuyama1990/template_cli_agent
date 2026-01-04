from abc import ABC

from mlip_autopipec.interfaces import ISampler


class BaseSampler(ISampler, ABC):
    """Abstract base class for all samplers."""
