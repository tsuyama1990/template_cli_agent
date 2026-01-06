from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

import py3Dmol
import streamlit as st
import yaml
from ase.db import connect
from stpy3mol import showmol


def run_pipeline_subprocess(config: dict):
    """Runs the main CLI pipeline in a subprocess."""
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".yaml") as tmp:
        yaml.dump(config, tmp)
        config_path = tmp.name

    command = ["mlip-autopipec", "run", "--config", config_path]
    process = subprocess.run(command, capture_output=True, text=True, check=False)
    st.text(process.stdout)
    if process.returncode != 0:
        st.error(process.stderr)
        return None

    db_path = f"{config['project_name']}.db"
    return db_path


st.set_page_config(layout="wide")
st.title("MLIP-AutoPipe Web UI")

# --- Sidebar for Configuration ---
st.sidebar.header("Configuration")
project_name = st.sidebar.text_input("Project Name", "web_ui_run")
elements_str = st.sidebar.text_input("Elements (comma-separated)", "Si")
temperature = st.sidebar.number_input("MD Temperature (K)", value=300.0, min_value=0.0, step=50.0)
md_steps = st.sidebar.number_input("MD Steps", value=100, min_value=10, step=10)

# --- Session State Initialization ---
if "last_run_result_path" not in st.session_state:
    st.session_state.last_run_result_path = None

# --- Main Action Button ---
if st.sidebar.button("Run Pipeline", type="primary"):
    elements = [e.strip() for e in elements_str.split(",")]
    config_data = {
        "project_name": project_name,
        "system": {
            "elements": elements,
            "lattice": "diamond" if "Si" in elements else "fcc",
            "num_structures": 1,
        },
        "exploration": {
            "temperature": temperature,
            "md_steps": md_steps,
            "mlip_model": "emt",
        },
        "sampling": {"method": "random", "fraction": 1.0},
    }

    with st.status("Running MLIP-AutoPipe pipeline...", expanded=True) as status:
        st.write("Executing CLI command...")
        db_path = run_pipeline_subprocess(config_data)
        if db_path:
            st.session_state.last_run_result_path = db_path
            status.update(label="Pipeline run completed!", state="complete", expanded=False)
        else:
            st.session_state.last_run_result_path = None
            status.update(label="Pipeline run failed!", state="error", expanded=True)

# --- Results Display ---
st.header("Results")
if st.session_state.last_run_result_path:
    db_path = Path(st.session_state.last_run_result_path)
    if db_path.exists():
        st.success(f"Output database created at: `{db_path}`")
        try:
            with connect(db_path) as db:
                if len(db) > 0:
                    atoms = db.get_atoms(id=1)
                    st.write("3D Visualization of a generated structure:")

                    xyz_view = py3Dmol.view(width=800, height=400)
                    xyz_view.addModel(
                        atoms.positions.tolist(),
                        "xyz",
                        {"vibrate": {"frames": 10, "amplitude": 0.1}},
                    )
                    xyz_view.setStyle({"stick": {}})
                    xyz_view.zoomTo()
                    showmol(xyz_view, height=400)
                else:
                    st.warning("Output database is empty.")
        except Exception as e:
            st.error(f"Failed to read or visualize results: {e}")
    else:
        st.error("Result database not found. The pipeline might have failed.")
else:
    st.info("Run the pipeline from the sidebar to see results.")
