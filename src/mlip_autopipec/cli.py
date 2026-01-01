import click
import yaml
from pydantic import ValidationError

from mlip_autopipec.config import load_config
from mlip_autopipec.workflow import run_cycle01_workflow


@click.group()
def cli():
    """MLIP-AutoPipe: Automated MLIP Generation."""
    pass


@cli.command(name="run-cycle")
@click.option(
    '--config',
    'config_path',
    required=True,
    type=click.Path(exists=True),
    help='Path to the YAML configuration file.'
)
def run_cycle(config_path: str):
    """
    Run a full labeling and training cycle.
    """
    try:
        config = load_config(config_path)
    except (FileNotFoundError, yaml.YAMLError, ValidationError) as e:
        click.echo(f"Error loading configuration: {e}", err=True)
        raise SystemExit(1) from e
    except Exception as e:
        click.echo(f"An unexpected error occurred: {e}", err=True)
        raise SystemExit(1) from e

    run_cycle01_workflow(config)


if __name__ == '__main__':
    cli()
