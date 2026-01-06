"""
Defines abstract interfaces for dependency injection.

This module contains abstract base classes that define the contracts for various
services and components in the application. Using these interfaces allows for
decoupling of components and easier testing via mock implementations.
"""

from .generators.base import BaseStructureGenerator
from .samplers.base import BaseSampler

# Re-exporting for a central location, as per the architecture diagram's intent.
# This makes it clear where to look for the system's contracts.
IGenerator = BaseStructureGenerator
ISampler = BaseSampler
