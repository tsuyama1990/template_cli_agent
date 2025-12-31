import click
import yaml

from .config.models import FullConfig
from .orchestrator import Orchestrator


@click.group()
def cli():
    """MLIP-AutoPipe command-line interface."""
    pass


@cli.command()
@click.option(
    "--config",
    "-c",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the input YAML configuration file.",
)
def run(config_path: str):
    """
    Runs the MLIP-AutoPipe workflow.
    """
    print(f"Loading configuration from: {config_path}")
    try:
        with open(config_path) as f:
            config_dict = yaml.safe_load(f)

        # Validate the configuration using Pydantic
        config = FullConfig(**config_dict)

    except yaml.YAMLError as e:
        click.echo(f"Error parsing YAML file: {e}", err=True)
        raise click.Abort() from e
    except Exception as e:
        # Pydantic's ValidationError will be caught here
        click.echo(f"Configuration validation error: {e}", err=True)
        raise click.Abort() from e

    # Instantiate and run the orchestrator
    orchestrator = Orchestrator(config=config)
    orchestrator.run_cycle_01()


if __name__ == "__main__":
    cli()
