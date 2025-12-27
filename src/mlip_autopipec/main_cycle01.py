import click
from mlip_autopipec.orchestrator_cycle01 import run_cycle01_workflow

@click.command()
@click.option('--config', 'config_path', type=click.Path(exists=True), required=True, help='Path to the YAML configuration file.')
@click.option('--structure', 'structure_path', type=click.Path(exists=True), required=True, help='Path to the atomic structure file (e.g., CIF, POSCAR).')
def main(config_path: str, structure_path: str):
    """
    Command-line interface for running the MLIP-AutoPipe Cycle 01 workflow.
    """
    run_cycle01_workflow(config_path, structure_path)

if __name__ == '__main__':
    main()
