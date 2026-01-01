from pathlib import Path

import click
import yaml
from mlip_autopipec.services.config_expander import ConfigExpander
from mlip_autopipec.workflow import initialize_and_run_workflow


@click.group()
def cli():
    """MLIP-AutoPipe: Automated MLIP Generation."""
    pass


@cli.command(name="run")
@click.option(
    "--input",
    "input_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the minimal `input.yaml` file.",
)
@click.option(
    "--output-config",
    "output_config_path",
    default="exec_config_dump.yaml",
    type=click.Path(dir_okay=False, writable=True),
    help="Path to write the full expanded configuration.",
)
def run(input_path: str, output_config_path: str):
    """
    Run the full `Generate -> Label -> Train` workflow from a minimal input.
    """
    try:
        # Define the path to the heuristic defaults file relative to the package
        defaults_path = Path(__file__).parent / "configs" / "sssp_defaults.json"

        # 1. Expand the minimal config into a full config
        click.echo("--- Stage 1: Expanding Configuration ---")
        expander = ConfigExpander(defaults_path=defaults_path)
        full_config = expander.expand_config(input_path=Path(input_path))

        # 2. Save the full config for user inspection and reproducibility
        with open(output_config_path, "w") as f:
            # Use Pydantic's model_dump to get a serializable dictionary
            config_dict = full_config.model_dump()
            yaml.dump(config_dict, f, sort_keys=False)
        click.echo(f"Full configuration saved to: {output_config_path}")

        # 3. Initialize and run the main workflow with the full config
        click.echo("\n--- Stage 2: Initializing and Running Main Workflow ---")
        initialize_and_run_workflow(config=full_config)

        click.echo("\nWorkflow finished successfully!", color="green")

    except (FileNotFoundError, yaml.YAMLError, ValueError) as e:
        click.echo(f"Error during workflow execution: {e}", err=True)
        raise SystemExit(1) from e
    except Exception as e:
        click.echo(f"An unexpected critical error occurred: {e}", err=True)
        raise SystemExit(1) from e


if __name__ == "__main__":
    cli()
