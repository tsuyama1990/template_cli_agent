from dataclasses import dataclass

from ac_cdd_core.presentation import ConsolePresenter
from ac_cdd_core.services.artifacts import ArtifactManager
from ac_cdd_core.services.contracts import ContractManager
from ac_cdd_core.services.file_ops import FilePatcher
from ac_cdd_core.sandbox import SandboxRunner


@dataclass
class ServiceContainer:
    file_patcher: FilePatcher
    contract_manager: ContractManager
    artifact_manager: ArtifactManager
    presenter: ConsolePresenter
    sandbox: SandboxRunner

    @classmethod
    def default(cls) -> "ServiceContainer":
        return cls(
            file_patcher=FilePatcher(),
            contract_manager=ContractManager(),
            artifact_manager=ArtifactManager(),
            presenter=ConsolePresenter(),
            sandbox=SandboxRunner(),
        )
