FROM python:3.12-slim

# Install system dependencies
# git and gh (GitHub CLI) are required
RUN apt-get update && apt-get install -y \
    git \
    curl \
    gnupg \
    && mkdir -p -m 755 /etc/apt/keyrings \
    && curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | gpg --dearmor -o /etc/apt/keyrings/githubcli-archive-keyring.gpg \
    && chmod go+r /etc/apt/keyrings/githubcli-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
    && apt-get update \
    && apt-get install -y gh \
    && rm -rf /var/lib/apt/lists/*

# Install uv
COPY --from=ghcr.io/astral-sh/uv:latest /uv /bin/uv

# Set working directory
WORKDIR /app

# Create directory for internal templates
RUN mkdir -p /opt/ac-cdd/templates

# Copy project files
# We assume the build context is the repository root
COPY pyproject.toml .
COPY dev_src/ ./dev_src/

# Install dependencies
# We install the package in editable mode or just dependencies
# Assuming we want to install the tool globally in the container
RUN uv pip install --system .

# Copy system prompts to internal template directory
# We assume they are in dev_documents/system_prompts based on previous findings
# Note: The prompt mentioned /opt/ac_cdd/templates (underscore vs hyphen).
# I used hyphen in config (/opt/ac-cdd/templates). I will stick to hyphen or check what I wrote.
# I wrote /opt/ac-cdd/templates in config.py default.
COPY dev_documents/system_prompts/*.md /opt/ac-cdd/templates/

# Entrypoint script
COPY entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/entrypoint.sh

# Environment variables
ENV PYTHONUNBUFFERED=1
ENV AC_CDD_TEMPLATE_PATH=/opt/ac-cdd/templates
# Ensure ac_cdd_config.py is looked for in CWD (/app)
# AC_CDD_CONFIG_PATH can be set by user if needed, but default is CWD.

ENTRYPOINT ["entrypoint.sh"]
CMD ["ac-cdd", "--help"]
