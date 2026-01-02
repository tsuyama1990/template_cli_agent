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

# Set working directory for the tool installation
# We install the tool into /opt/ac_cdd/ac_cdd_core
# Although previously we copied to /opt/ac_cdd/ac_cdd_core and installed .,
# now the package root is ac_cdd_core (and src for user).
# The build context root has pyproject.toml and ac_cdd_core/.

WORKDIR /opt/ac_cdd

# Copy project files for the tool
COPY pyproject.toml .
# Copy the core package directory
COPY ac_cdd_core/ ./ac_cdd_core/

# Install the tool and its dependencies into the system python environment
# Note: we are installing from . which contains pyproject.toml
RUN uv pip install --system .

# Copy system prompts to internal template directory
RUN mkdir -p /opt/ac_cdd/templates
COPY dev_documents/system_prompts/*.md /opt/ac_cdd/templates/

# Create /app for user mount
WORKDIR /app

# Copy entrypoint
COPY entrypoint.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/entrypoint.sh

# Environment variables
ENV PYTHONUNBUFFERED=1
ENV AC_CDD_TEMPLATE_PATH=/opt/ac_cdd/templates

ENTRYPOINT ["ac-cdd"]
CMD ["--help"]
