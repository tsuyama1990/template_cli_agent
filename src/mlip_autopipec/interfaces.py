from abc import ABC, abstractmethod

from ase import Atoms

from mlip_autopipec.config import DFTResult


class ILabelingEngine(ABC):
    """Interface for a labeling engine that processes atomic structures."""

    @abstractmethod
    def label_structure(self, atoms: Atoms) -> DFTResult:
        """
        Labels a single atomic structure by running a calculation.

        Args:
            atoms: The ASE Atoms object to label.

        Returns:
            A DFTResult object containing the calculated energy, forces, and stress.
        """
        pass


class ITrainingEngine(ABC):
    """Interface for a training engine that trains an MLIP model."""

    @abstractmethod
    def train(self, training_data: list[tuple[Atoms, DFTResult]]) -> None:
        """
        Trains the MLIP model using the provided labeled data.

        Args:
            training_data: A list of tuples, each containing an Atoms object
                           and its corresponding DFTResult.
        """
        pass


class IStructureGenerator(ABC):
    """Interface for a structure generator."""

    @abstractmethod
    def generate(self) -> list[Atoms]:
        """
        Generates a list of initial atomic structures.

        Returns:
            A list of ase.Atoms objects.
        """
        pass


class IWorkflowOrchestrator(ABC):
    """Interface for the main workflow orchestrator."""

    @abstractmethod
    def run(self) -> None:
        """Runs the main workflow."""
        pass

    @abstractmethod
    def label_structure_by_id(self, structure_id: int) -> None:
        """
        Runs the labeling process for a specific structure ID.

        Args:
            structure_id: The ID of the structure to label.
        """
        pass

    @abstractmethod
    def run_training(self) -> None:
        """Runs the MLIP model training process."""
        pass


class IProcessRunner(ABC):
    """Interface for running external processes."""

    @abstractmethod
    def run(self, command: list[str], stdout_path: str) -> None:
        """
        Runs an external command and redirects its stdout to a file.

        Args:
            command: The command to run as a list of strings.
            stdout_path: The path to the file to redirect stdout to.

        Raises:
            subprocess.CalledProcessError: If the command returns a non-zero exit code.
        """
        pass
