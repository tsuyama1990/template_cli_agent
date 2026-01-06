"""
Custom exceptions for the MLIP-AutoPipe application.

Defining custom exceptions allows for more specific error handling and clearer
error messages, which improves the robustness and user-friendliness of the tool.
"""


class PipelineError(Exception):
    """Base class for all custom exceptions in this application."""



class ConfigurationError(PipelineError):
    """Raised for errors in the user's configuration file."""



class GenerationError(PipelineError):
    """Raised for errors during the structure generation stage."""



class ExplorationError(PipelineError):
    """Raised for errors during the MD exploration stage."""



class SamplingError(PipelineError):
    """Raised for errors during the sampling stage."""



class StorageError(PipelineError):
    """Raised for errors during the database storage stage."""
