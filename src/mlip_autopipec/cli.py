from pathlib import Path

import click
import yaml

from mlip_autopipec.data.models import Cycle01Config
from mlip_autopipec.orchestrator import Orchestrator


@click.group()
def cli():
    """MLIP-AutoPipe command-line interface."""
    pass


@cli.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Path to the cycle configuration YAML file.",
)
def run_cycle(config_path: Path):
    """Run a full MLIP-AutoPipe cycle."""
    with open(config_path) as f:
        config_dict = yaml.safe_load(f)

    # Check if the database file exists, create if not
    db_path_str = config_dict.get("database_path", "default.db")
    db_path = Path(db_path_str)
    if not db_path.exists():
        print(f"Database not found at {db_path}, creating a new one.")
        db_path.touch()

    config = Cycle01Config(**config_dict)

    orchestrator = Orchestrator(config)
    orchestrator.run_label_and_train_workflow()

    click.echo("Cycle finished successfully. Trained model saved to trained_model.pt")


if __name__ == "__main__":
    cli()
