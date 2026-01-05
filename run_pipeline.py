import yaml
from mlip_autopipec.config.models import FullConfig
from mlip_autopipec.core.orchestrator import PipelineRunner

def main():
    config_path = "config_valid.yaml"
    with open(config_path) as f:
        config_dict = yaml.safe_load(f)

    full_config = FullConfig(**config_dict)
    runner = PipelineRunner(full_config)
    runner.run()

if __name__ == "__main__":
    main()
