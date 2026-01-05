"""Unit tests for the configuration models."""
import yaml
import pytest
from pydantic import ValidationError
from mlip_autopipec.config.models import FullConfig

VALID_CONFIG_YAML = """
project_name: test_project
system:
  elements: ['Fe', 'Pt']
  composition: {'Fe': 0.5, 'Pt': 0.5}
  lattice: 'fcc'
  num_structures: 10
exploration:
  temperature: 300.0
sampling:
  method: 'random'
  fraction: 0.8
"""

def test_valid_config():
    """Test that a valid configuration is parsed correctly."""
    config_dict = yaml.safe_load(VALID_CONFIG_YAML)
    config = FullConfig(**config_dict)

    assert config.project_name == "test_project"
    assert config.system.elements == ['Fe', 'Pt']
    assert config.system.composition == {'Fe': 0.5, 'Pt': 0.5}
    assert config.system.lattice == 'fcc'
    assert config.system.num_structures == 10
    assert config.exploration.temperature == 300.0
    assert config.sampling.method == 'random'
    assert config.sampling.fraction == 0.8

def test_invalid_composition_sum():
    """Test that a composition not summing to 1.0 raises an error."""
    config_dict = yaml.safe_load(VALID_CONFIG_YAML)
    config_dict['system']['composition'] = {'Fe': 0.6, 'Pt': 0.5}
    with pytest.raises(ValidationError, match="Composition probabilities must sum to 1.0."):
        FullConfig(**config_dict)

def test_mismatched_elements():
    """Test that mismatched elements and composition keys raise an error."""
    config_dict = yaml.safe_load(VALID_CONFIG_YAML)
    config_dict['system']['composition'] = {'Fe': 0.5, 'Au': 0.5}
    with pytest.raises(ValidationError, match="Elements and composition keys must match."):
        FullConfig(**config_dict)

def test_negative_temperature():
    """Test that a negative temperature raises an error."""
    config_dict = yaml.safe_load(VALID_CONFIG_YAML)
    config_dict['exploration']['temperature'] = -100.0
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        FullConfig(**config_dict)

def test_invalid_sampling_fraction():
    """Test that a sampling fraction > 1.0 raises an error."""
    config_dict = yaml.safe_load(VALID_CONFIG_YAML)
    config_dict['sampling']['fraction'] = 1.1
    with pytest.raises(ValidationError, match="Input should be less than or equal to 1"):
        FullConfig(**config_dict)

def test_incorrect_data_type():
    """Test that an incorrect data type raises an error."""
    config_dict = yaml.safe_load(VALID_CONFIG_YAML)
    config_dict['system']['num_structures'] = "ten"
    with pytest.raises(ValidationError, match="Input should be a valid integer"):
        FullConfig(**config_dict)
