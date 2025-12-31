import click
from ase import Atoms
from mlip_autopipec.orchestrator import Orchestrator

@click.group()
def cli():
    """MLIP-AutoPipe Command Line Interface"""
    pass

@cli.command()
@click.option("--qe-command", default="pw.x", help="Command to run Quantum Espresso.")
@click.option("--pseudo-dir", default="./pseudos", help="Directory with pseudopotentials.")
@click.option("--db-path", default="mlip.db", help="Path to the ASE database.")
def run_cycle01(qe_command, pseudo_dir, db_path):
    """
    Runs the simple, hardcoded workflow for Cycle 01.
    """
    # For Cycle 1, we use a hardcoded H2 molecule as the input.
    atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]])

    orchestrator = Orchestrator(
        db_path=db_path,
        qe_command=qe_command,
        pseudo_dir=pseudo_dir
    )
    orchestrator.run_cycle01_workflow(atoms)
    click.echo("Cycle 01 workflow finished.")

if __name__ == "__main__":
    cli()
