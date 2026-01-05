# Cycle 5 Specification: Web UI and User Experience

## 1. Summary

Cycle 5 is the final planned cycle for the MLIP-AutoPipe project and focuses entirely on user experience and accessibility. While the powerful CLI-driven pipeline is the core of the application, a graphical user interface is essential for lowering the barrier to entry, facilitating interactive exploration, and appealing to a broader user base. This cycle will deliver a web-based UI that wraps the core pipeline, providing a user-friendly way to configure, execute, and monitor the data generation process.

The key deliverables for this cycle are:
1.  **A Web-Based User Interface:** A web application will be developed, likely using a Python framework like Streamlit or FastAPI with a simple frontend, that serves as a graphical front-end to the `mlip_autopipec` library.
2.  **Interactive Configuration Builder:** Instead of manually writing YAML files, users will be guided through a series of interactive widgets (sliders, dropdowns, text boxes) to build a valid configuration for the pipeline. This reduces errors and makes the configuration process more intuitive.
3.  **Pipeline Monitoring:** The UI will provide real-time (or near-real-time) feedback on the status of a running pipeline. This includes showing which stage is currently active, displaying log messages, and indicating overall progress.
4.  **Results Visualization:** The UI will include a component for visualizing the outputs. This will feature an interactive 3D viewer to inspect the atomic structures generated and simple plots to show key data, such as the distribution of potential energies in the final dataset.

By the end of Cycle 5, the MLIP-AutoPipe project will be a complete and polished application. It will offer two distinct interfaces—a powerful CLI for experts and automation, and an intuitive Web UI for new users and interactive work. This dual-interface approach maximizes the project's utility and impact, making it accessible to everyone from computational power-users to experimental scientists who are new to simulation.

## 2. System Architecture

This cycle introduces a new entry point to the application, `main_gui.py`, and the associated UI code. It will reuse all of the backend library code developed in Cycles 1-4.

```
src/mlip_autopipec/
├── __init__.py
├── cli.py
├── **config.py**               # Pydantic models may be enhanced with UI-specific metadata
├── **main_gui.py**             # **New entry point** for the Web UI
├── database/
│   └── ase_db_wrapper.py
├── domain/
│   └── ...
├── engines/
│   └── ...
├── explorers/
│   └── ...
├── generators/
│   └── ...
├── samplers/
│   └── ...
├── **ui/**
│   ├── **__init__.py**
│   ├── **app.py**                # Main application logic for the UI
│   ├── **components.py**         # Reusable UI components (e.g., 3D viewer)
│   └── **state.py**              # UI state management
└── utils/
│   └── ...
└── workflow_orchestrator.py
```

**File Blueprints:**

*   **`main_gui.py`**:
    *   A simple script that imports the main application object from `ui/app.py` and launches it. For Streamlit, this might just be a few lines that import and run the main script. For FastAPI, it would launch the Uvicorn server.

*   **`ui/app.py`**:
    *   This will be the core of the web application.
    *   It will be responsible for laying out the main user interface, including sections for "Configuration", "Execution", and "Results".
    *   **Configuration Section:** It will use UI widgets to build up a Python dictionary or a Pydantic `MainConfig` object in memory. For example, a dropdown to select the generator mode (`alloy`/`ionic`), sliders for temperature, and text inputs for elements.
    *   **Execution Section:** A "Run Pipeline" button will trigger the main workflow. When clicked, it will take the configuration object, instantiate the `WorkflowOrchestrator` from the existing library code, and run it in a separate thread or process to avoid blocking the UI.
    *   **Monitoring:** It will display the standard output/logs from the running pipeline process in a text area on the page. It will also show a progress bar to indicate the overall status.

*   **`ui/components.py`**:
    *   This file will contain reusable UI modules to keep the main `app.py` file clean.
    *   A `structure_viewer(atoms)` function will be created that takes an `ase.Atoms` object and renders it using an embedded 3D viewer library like `py3Dmol`.
    *   A `plot_energy_distribution(db_path)` function that connects to a database, extracts the energies, and uses a library like Matplotlib or Plotly to generate and display a histogram.

*   **`ui/state.py`**:
    *   Web applications are stateless, so this file will manage the session state. For Streamlit, this involves using the `st.session_state` object.
    *   It will be responsible for storing the current configuration, the status of the pipeline (e.g., "idle", "running", "finished"), and the path to the output database.

## 3. Design Architecture

The design of the UI prioritizes a clean separation between the user interface (the "view") and the backend application logic (the "model").

**Design Principles:**

*   **Reusability of Core Logic:** The most important principle is that the UI is just another "client" of the existing `mlip_autopipec` library. It should **not** reimplement any of the core pipeline logic. It will import and use the `WorkflowOrchestrator`, `AseDBWrapper`, and other components directly. This ensures that both the CLI and the Web UI are always running the exact same underlying code, which dramatically reduces the maintenance burden.
*   **Non-Blocking Execution:** Running the pipeline is a long, computationally intensive task. It is unacceptable for the UI to freeze while the backend is running. The design must therefore run the `WorkflowOrchestrator` in a background process or thread. The UI will then periodically poll for status updates (e.g., by reading from a log file or a shared queue) to update its monitoring components. This asynchronous design is essential for a responsive user experience.
*   **Component-Based UI:** The UI will be built from small, reusable components (defined in `ui/components.py`). For example, the 3D structure viewer is a component. This makes the code more modular and easier to test and debug. If we decide to switch from `py3Dmol` to a different viewer library, we only need to change one function in one file.
*   **State Management:** The UI state (like the user's current configuration choices) will be explicitly managed in `ui/state.py`. This centralizes state-related logic and makes it easier to understand how data flows through the application as the user interacts with it.

## 4. Implementation Approach

1.  **Framework Selection and Setup:**
    *   Choose a UI framework. Streamlit is an excellent choice for its speed of development and simplicity, making it ideal for this kind of data-centric application.
    *   Add the chosen framework (`streamlit`) and any required visualization libraries (`py3dmol`, `nglview`, `plotly`) as dependencies in `pyproject.toml`.

2.  **Configuration UI (`ui/app.py`):**
    *   Create the main application layout.
    *   Use `st.sidebar` to create a configuration panel.
    *   Add widgets for all the main configuration options defined in the Pydantic models in `config.py`. Use `st.selectbox` for choices like generator mode, `st.slider` for numerical values like temperature, and `st.text_input` for element lists.
    *   As the user changes the widgets, update the application's session state.

3.  **Backend Integration (`ui/app.py`):**
    *   Create a "Run Pipeline" button.
    *   Write the `on_click` handler for this button. This function will:
        a. Retrieve the configuration from the session state.
        b. Convert it into the Pydantic `MainConfig` object.
        c. Use Python's `multiprocessing` or `subprocess` module to launch a new process that runs the `WorkflowOrchestrator`. The logs from this process should be redirected to a file.
        d. Store the process handle in the session state and set the application status to "running".

4.  **Monitoring UI (`ui/app.py`):**
    *   Create a placeholder on the main page for log output.
    *   Implement a mechanism that, while the status is "running", periodically (e.g., every 2 seconds) reads the latest content from the log file and displays it in the placeholder. This gives the user real-time feedback.
    *   Add a `st.progress` bar that is updated as the pipeline moves from one stage to the next.

5.  **Results Visualization (`ui/components.py`, `ui/app.py`):**
    *   Create the `plot_energy_distribution` component. It will take the path to the output database, connect to it, fetch the data, and render a plot.
    *   Create the `structure_viewer` component.
    *   Once the pipeline status is "finished", the main app will display a "Results" section. It will use the `AseDBWrapper` to fetch a few sample structures from the final dataset and display them using the `structure_viewer`. It will also call the `plot_energy_distribution` component to show the summary graphic.

## 5. Test Strategy

Testing a UI is different from testing backend logic. The focus is on the user's interaction flow and the visual presentation of data.

**Unit Testing Approach (Min 300 words):**

Unit testing UIs can be complex. We will focus on testing the backend logic of the UI components, not their graphical rendering.
*   **Component Logic:** We can test the functions in `ui/components.py` in isolation.
    *   For `plot_energy_distribution`, we can create a dummy database file with known energy values. We will call the function with the path to this dummy file. Instead of asserting the plot image, we can have the function return the data it prepared for plotting and assert that this data is correct.
    *   For `structure_viewer`, we can test that it correctly converts an `ase.Atoms` object into the specific string or dictionary format expected by the `py3Dmol` library.
*   **State Management:** We can test the state transition logic. For example, we can write a test that checks if clicking the "Run" button correctly changes the application state in `st.session_state` from "idle" to "running". This can be done by calling the callback function directly in the test and then inspecting a mock `session_state` object.

**Integration Testing Approach (Min 300 words):**

End-to-end testing is the most effective way to validate the UI. We will use a browser automation framework to simulate a real user.
*   **Tool Selection:** We will use a tool like Playwright, which allows us to write Python scripts that control a web browser.
*   **Test Scenario (Full Workflow):**
    1.  **Setup:** The test script will first launch the Streamlit application using `subprocess`.
    2.  **Interaction:** The Playwright script will then:
        a. Connect to the running application in a browser.
        b. Interact with the widgets to fill out the configuration form (e.g., `page.select_option`, `page.fill`, `page.click`).
        c. Click the "Run Pipeline" button.
    3.  **Monitoring:** The script will wait and assert that the UI updates to show a "running" status and that the log text area starts to fill with content.
    4.  **Verification:** The script will wait until a "Finished" message appears in the UI. It will then:
        a. Assert that a plot canvas is now visible on the page.
        b. Assert that a 3D viewer element is now present.
        c. It can even take a screenshot of the final page (`page.screenshot`) and save it as a test artifact.
    This automated, end-to-end test provides very high confidence that the entire UI workflow is functioning correctly, from user input to the final display of results. It will be the primary UAT for this cycle.
