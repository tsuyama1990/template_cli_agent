"""
This module defines the command-line interface (CLI) for the application.

It uses the `click` library to create a user-friendly interface for executing
the main workflows, such as labeling structures and training models. The CLI
is responsible for parsing user input, orchestrating the configuration
expansion, and invoking the core application logic.
"""

import click

from mlip_autopipec.app import Application
from mlip_autopipec.config import ConfigExpander, UserInputConfig


@click.group()
def app() -> None:
    """MLIP-AutoPipe: Automated MLIP Generation."""
    pass


@app.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the user input YAML file.",
)
@click.option(
    "--db-path",
    type=str,
    default=None,
    help="Override the path to the ASE database file specified in the config.",
)
def run(config_path: str, db_path: str | None) -> None:
    """Runs the full MLIP generation workflow."""
    user_input = UserInputConfig.from_yaml(config_path)
    expander = ConfigExpander()
    full_config = expander.expand(user_input)

    if db_path:
        full_config.db_path = db_path

    full_config.to_yaml("exec_config_dump.yaml")

    application = Application(full_config)
    application.run()


@app.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the user input YAML file.",
)
@click.option(
    "--id",
    required=True,
    type=int,
    help="The ID of the structure to label.",
)
@click.option(
    "--db-path",
    type=str,
    default=None,
    help="Override the path to the ASE database file specified in the config.",
)
def label(config_path: str, id: int, db_path: str | None) -> None:
    """Runs the DFT labeling for a single structure."""
    user_input = UserInputConfig.from_yaml(config_path)
    expander = ConfigExpander()
    full_config = expander.expand(user_input)

    if db_path:
        full_config.db_path = db_path

    full_config.to_yaml("exec_config_dump.yaml")

    application = Application(full_config)
    application.label_structure(id)


@app.command()
@click.option(
    "--config",
    "config_path",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
    help="Path to the user input YAML file.",
)
@click.option(
    "--db-path",
    type=str,
    default=None,
    help="Override the path to the ASE database file specified in the config.",
)
def train(config_path: str, db_path: str | None) -> None:
    """Trains the MLIP model on the existing labeled data."""
    user_input = UserInputConfig.from_yaml(config_path)
    expander = ConfigExpander()
    full_config = expander.expand(user_input)

    if db_path:
        full_config.db_path = db_path

    full_config.to_yaml("exec_config_dump.yaml")

    application = Application(full_config)
    application.train()


if __name__ == "__main__":
    app()
