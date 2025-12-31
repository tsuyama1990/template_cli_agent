import click
import yaml
from pydantic import ValidationError
import os

from mlip_autopipec.orchestrator import Orchestrator
from mlip_autopipec.data.models import Cycle01Config

@click.group()
def cli():
    """MLIP-AutoPipe: Automated MLIP Generation."""
    pass

@cli.command(name="run-cycle")
@click.option('--config', 'config_path', required=True, type=click.Path(exists=True), help='Path to the YAML configuration file.')
def run_cycle(config_path: str):
    """
    Run a full labeling and training cycle.
    """
    try:
        print(f"Loading configuration from: {config_path}")
        with open(config_path, 'r') as f:
            config_dict = yaml.safe_load(f)

        # Resolve the database path relative to the config file location
        config_dir = os.path.dirname(os.path.abspath(config_path))
        db_path = os.path.join(config_dir, config_dict.get('database_path', ''))
        config_dict['database_path'] = db_path

        config = Cycle01Config(**config_dict)

    except FileNotFoundError:
        click.echo(f"Error: Configuration file not found at '{config_path}'", err=True)
        raise SystemExit(1)
    except yaml.YAMLError as e:
        click.echo(f"Error: Could not parse YAML file: {e}", err=True)
        raise SystemExit(1)
    except ValidationError as e:
        click.echo(f"Error: Configuration validation failed:\n{e}", err=True)
        raise SystemExit(1)
    except Exception as e:
        click.echo(f"An unexpected error occurred: {e}", err=True)
        raise SystemExit(1)

    orchestrator = Orchestrator(config)
    orchestrator.run_label_and_train_workflow()

if __name__ == '__main__':
    cli()
