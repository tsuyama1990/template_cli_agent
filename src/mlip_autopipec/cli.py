# src/mlip_autopipec/cli.py

import click
import yaml
from pydantic import ValidationError

from mlip_autopipec.configs.models import MainConfig
from mlip_autopipec.workflow import WorkflowOrchestrator


@click.command()
@click.option(
    "--config-file",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the main YAML configuration file.",
)
@click.option(
    "--database-file",
    required=True,
    type=click.Path(dir_okay=False),
    help="Path to the ASE database file to use or create.",
)
@click.option(
    "--input-file",
    type=click.Path(exists=True, dir_okay=False),
    help="Optional path to an XYZ file with initial structures to add to the DB.",
)
def main(config_file: str, database_file: str, input_file: str):
    """
    Main entry point for the MLIP-AutoPipe CLI.
    """
    try:
        with open(config_file) as f:
            raw_config = yaml.safe_load(f)
        config = MainConfig(**raw_config)
    except FileNotFoundError:
        click.echo(f"Error: Configuration file not found at {config_file}", err=True)
        return
    except yaml.YAMLError as e:
        click.echo(f"Error parsing YAML configuration: {e}", err=True)
        return
    except ValidationError as e:
        click.echo(f"Configuration validation error:\n{e}", err=True)
        return

    orchestrator = WorkflowOrchestrator(config)
    orchestrator.run(database_path=database_file, input_file_path=input_file)


if __name__ == "__main__":
    main()
