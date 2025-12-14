import json
from typing import Any

from google import genai
from google.genai import types
from rich.console import Console

from ac_cdd.agent_interface import AgentInterface
from ac_cdd.utils import logger

# Constants for cost estimation
# Placeholder rates: Input $0.10/1M, Output $0.40/1M
COST_PER_1M_INPUT_TOKENS = 0.10
COST_PER_1M_OUTPUT_TOKENS = 0.40

console = Console()

class GeminiApiClient(AgentInterface):
    """
    Agent interface implementation for Gemini via google-genai library.
    """
    def __init__(self, api_key: str, model: str = "gemini-pro") -> None:
        self.api_key = api_key
        self.model = model
        self.client = genai.Client(api_key=self.api_key)

    def _calculate_cost(self, usage_metadata: Any) -> float:
        """
        Calculates the estimated cost based on token usage.
        usage_metadata object from google-genai
        """
        if not usage_metadata:
            return 0.0

        # Access attributes directly as per google-genai types
        prompt_tokens = getattr(usage_metadata, "prompt_token_count", 0)
        candidates_tokens = getattr(usage_metadata, "candidates_token_count", 0)

        input_cost = (prompt_tokens / 1_000_000) * COST_PER_1M_INPUT_TOKENS
        output_cost = (candidates_tokens / 1_000_000) * COST_PER_1M_OUTPUT_TOKENS

        return input_cost + output_cost

    def _log_cost(self, usage_metadata: Any) -> None:
        """
        Logs the estimated cost using Rich.
        """
        cost = self._calculate_cost(usage_metadata)

        prompt_tokens = getattr(usage_metadata, "prompt_token_count", 0)
        candidates_tokens = getattr(usage_metadata, "candidates_token_count", 0)

        console.print(f"[yellow]💰 Est. Cost: ${cost:.6f}[/yellow]")
        logger.info(f"Token Usage - Prompt: {prompt_tokens}, Candidates: {candidates_tokens}")

    def start_task(self, prompt: str, **kwargs: Any) -> str:
        """
        Starts a new task (One-off generation for Gemini usually).
        If 'json_mode' is True in kwargs, enforces JSON response.
        """
        return self.generate_content(prompt, **kwargs)

    def send_message(self, prompt: str, **kwargs: Any) -> str:
        """
        Sends a message. For Gemini one-off, alias to start_task.
        """
        return self.generate_content(prompt, **kwargs)

    def generate_content(self, prompt: str, **kwargs: Any) -> str:
        json_mode = kwargs.get("json_mode", False)

        config_args: dict[str, Any] = {}
        if json_mode:
            config_args["response_mime_type"] = "application/json"

        # Add system instruction if provided in prompt or kwargs?
        # The user instruction says to prepend system prompt,
        # so we assume 'prompt' contains everything.

        try:
            response = self.client.models.generate_content(
                model=self.model,
                contents=prompt,
                config=types.GenerateContentConfig(**config_args)
            )

            if response.usage_metadata:
                self._log_cost(response.usage_metadata)

            return str(response.text)
        except Exception as e:
            logger.error(f"Gemini API Error: {e}")
            # Return JSON error structure if in JSON mode to avoid crashing parsers
            if json_mode:
                return json.dumps({"approved": False, "comments": [f"API Error: {str(e)}"]})
            raise

    def get_status(self) -> str:
        return "Gemini Client Ready"
