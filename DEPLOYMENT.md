# Deployment and Maintenance Strategy: MLIP-AutoPipe

This document outlines the strategy for deploying the MLIP-AutoPipe application and maintaining it post-deployment.

## 1. Deployment Strategy

The application consists of two main components with distinct deployment strategies: the core CLI tool and the Web UI.

### 1.1. CLI Tool Deployment

The primary user-facing component is the CLI tool. It will be packaged and distributed as a standard Python package on the Python Package Index (PyPI).

*   **Packaging:** We will use `hatchling`, as configured in `pyproject.toml`, to build the distributable wheel. The `project.scripts` section will define the `mlip-autopipec` command-line entrypoint.
*   **Distribution:**
    1.  Upon a new release (e.g., after a successful cycle), a version tag (e.g., `v0.1.0`) will be created in the Git repository.
    2.  A GitHub Actions workflow will be triggered by this tag.
    3.  The workflow will automatically build the Python wheel.
    4.  It will then publish the wheel to PyPI, making it installable via `pip install mlip-autopipec`.
*   **Target Audience:** This deployment method is suitable for researchers and developers who are comfortable with Python and command-line environments.

### 1.2. Web UI Deployment

The Web UI is a Streamlit application designed for broader accessibility. It will be deployed as a Docker container to ensure all dependencies are encapsulated and the environment is reproducible.

*   **Containerization:** A `Dockerfile` will be created specifically for the Web UI. This will:
    *   Start from a standard Python base image.
    *   Install all necessary Python packages, including `mlip-autopipec` itself.
    *   Set the `ENTRYPOINT` to run the Streamlit application (`streamlit run src/mlip_autopipec/main_gui.py`).
*   **Distribution:**
    1.  The Docker image will be automatically built and pushed to a container registry (e.g., Docker Hub, GitHub Container Registry) by the same GitHub Actions workflow used for the PyPI release.
    2.  The image will be tagged with the same version number.
*   **Deployment Scenarios:**
    *   **Local Use:** Users can easily run the UI on their local machine with a single `docker run` command.
    *   **Server Deployment:** The container can be deployed to any cloud provider (AWS, GCP, etc.) or on-premise server that supports Docker, allowing a research group to host a shared instance.

## 2. Post-Deployment Maintenance

A plan for maintenance is critical to ensure the long-term health and viability of the deployed application.

### 2.1. Monitoring

*   **Error Tracking:** For the Web UI, we will integrate a simple error tracking service (like Sentry's free tier) to capture and report any unhandled exceptions that occur in the deployed environment. This will provide immediate visibility into bugs affecting users.
*   **Usage Analytics (Optional):** If user adoption grows, we can consider adding opt-in, privacy-respecting analytics to understand which features are most used, helping to guide future development.

### 2.2. Updates and Patching

*   **Release Cadence:** New features will be released at the end of each successful development cycle.
*   **Bug Fixes:** Critical bugs will be addressed via hotfix releases. A hotfix will consist of a small, targeted change cherry-picked onto a release branch. It will go through the full testing and code review process before being deployed.
*   **Dependency Updates:** We will schedule a regular review (e.g., quarterly) of our external dependencies for any security vulnerabilities. Patches will be applied and tested before being rolled out in a new release.

### 2.3. Rollback Procedures

In the event that a deployment introduces a critical regression, we must be able to quickly revert to a stable state.

*   **CLI Tool:** Since the CLI is distributed via PyPI, users can roll back by simply installing a previous version: `pip install mlip-autopipec==<previous_version>`. Release notes will clearly document any major changes.
*   **Web UI:** The use of versioned Docker images makes rollbacks straightforward. If the `v0.2.0` image is found to be faulty, the deployment can be reverted to the `v0.1.0` image simply by stopping the new container and starting the old one. The containerized approach ensures this is a low-risk operation.
