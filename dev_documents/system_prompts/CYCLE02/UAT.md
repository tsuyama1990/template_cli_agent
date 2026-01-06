# UAT: CYCLE 02 - Advanced Features and UI

This document outlines the User Acceptance Testing (UAT) plan for the deliverables of Cycle 02. The goal is to verify the new scientific features (Ionic Generator, Hybrid MD/MC, ZBL, FPS) and the proof-of-concept Web UI.

This UAT will be performed using a Jupyter Notebook (`UAT_Cycle02.ipynb`) for scientific validation and manual testing for the Web UI.

## 1. Test Scenarios (Scientific Validation)

| Scenario ID | Test Scenario | Priority |
|---|---|---|
| UAT-C02-001 | **Ionic Structure Generation (Scientific Validation)** | High |
| UAT-C02-002 | **Verification of Hybrid MD/MC Exploration** | High |
| UAT-C02-003 | **Verification of ZBL Potential Stability** | High |
| UAT-C02-004 | **Verification of Farthest Point Sampling (FPS)** | Medium |

---

### **Scenario UAT-C02-001: Ionic Structure Generation (Scientific Validation)**

*   **Summary:** This scenario tests the `IonicGenerator`. The user will run the pipeline with a configuration for a simple ionic compound (e.g., NaCl). The notebook will then analyze the output database.
*   **Expected Outcome:** The pipeline completes successfully. The notebook analysis must show that **every structure in the database is charge-neutral**. It will also verify that the composition is correct (50% Na, 50% Cl).

---

### **Scenario UAT-C02-002: Verification of Hybrid MD/MC Exploration**

*   **Summary:** This scenario verifies that the Monte Carlo atom swaps are functioning correctly. The user will run a short exploration for an alloy system with MC swaps enabled.
*   **Expected Outcome:** The Jupyter notebook will load the final database and the initial seed structure. It will then assert that, while the final structures have different atomic arrangements from the initial one, the **overall chemical composition of every structure in the database remains identical** to the initial composition. This proves the MC swaps are correctly preserving the number of atoms of each element.

---

### **Scenario UAT-C02-003: Verification of ZBL Potential Stability**

*   **Summary:** This scenario demonstrates the stabilizing effect of the ZBL potential. The user will run two high-temperature simulations on the same initial structure: one with ZBL mixing disabled, and one with it enabled.
*   **Expected Outcome:** The run without ZBL should either fail with a physics violation error or result in a database containing structures with unrealistically high potential energies. The run **with ZBL enabled must complete successfully**, and the resulting database should contain structures with a distribution of energies that is physically plausible, which will be visualized in the notebook.

---

### **Scenario UAT-C02-004: Verification of Farthest Point Sampling (FPS)**

*   **Summary:** This scenario compares the output of FPS with the Random Sampler from Cycle 01. The user will run a long exploration to generate a large trajectory, and then sample from it using both methods.
*   **Expected Outcome:** The Jupyter notebook will calculate a diversity metric (e.g., the average pairwise distance in SOAP descriptor space) for the datasets produced by both Random Sampling and FPS. The notebook will then assert that the **diversity metric for the FPS-sampled set is significantly higher** than for the randomly sampled set, demonstrating the effectiveness of the advanced sampling algorithm.

## 2. Test Scenarios (Web UI)

| Scenario ID | Test Scenario | Priority |
|---|---|---|
| UAT-C02-005 | **Basic Web UI Functionality** | Medium |

---

### **Scenario UAT-C02-005: Basic Web UI Functionality**

*   **Summary:** This is a manual test of the proof-of-concept Web UI. The user will launch the Streamlit application.
*   **Expected Outcome:**
    1.  The user can successfully launch the Web UI.
    2.  The UI provides a mechanism to upload a `config.yaml` file.
    3.  The user can click a "Run Pipeline" button to start a data generation job as a background process.
    4.  The UI provides some visual feedback that the job is running.
    5.  After the job is complete, the user can view at least one of the generated structures in a 3D visualizer within the UI.
