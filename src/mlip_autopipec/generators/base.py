from abc import ABC

from mlip_autopipec.interfaces import IStructureGenerator


class BaseStructureGenerator(IStructureGenerator, ABC):
    """Abstract base class for all structure generators."""
