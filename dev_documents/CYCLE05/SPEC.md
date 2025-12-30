# CYCLE 05: SPEC.md

## 1. Summary

This document provides the detailed technical specifications for Cycle 05 of the Autonomous Development Environment (AC-CDD) project. This is the final and most advanced cycle, focusing on the implementation of the iterative fixing loop and the full orchestration of the end-to-end development workflow using LangGraph. This cycle will bring together all the components developed in the previous cycles to create a truly autonomous, self-correcting development pipeline. The core feature of this cycle is the ability for the system to not only detect issues, via the Auditor and QA Analyst agents, but to also autonomously fix them by re-engaging the Coder agent with the feedback from the quality assurance stage.

The implementation will involve creating a sophisticated orchestration graph using LangGraph, which will manage the state of the development cycle and control the flow of execution between the Coder, Auditor, and QA Analyst agents. This graph will be designed to handle the conditional logic required for the iterative fixing loop, allowing the system to cycle between coding, auditing, and testing until the code is of a high enough quality to be accepted.

The successful completion of this cycle will result in a fully realised version of the AC-CDD system, capable of taking high-level requirements and autonomously producing high-quality, fully tested code with minimal human intervention. This will represent the culmination of all the work done in the project, delivering on the promise of a truly autonomous development environment.

## 2. System Architecture

The system architecture for Cycle 05 is centred around the introduction of a powerful orchestration engine built with LangGraph.

The new and updated components are:

*   **Orchestration Engine (LangGraph):** The core of this cycle is the replacement of the simple, linear orchestration logic in the `CoderService` with a sophisticated, stateful orchestration graph built using LangGraph. This graph will manage the entire `run-cycle` workflow, including the iterative fixing loop.

*   **Coder Service:** The `CoderService` will be refactored to use the new LangGraph-based Orchestration Engine. It will be responsible for initiating the graph and providing it with the necessary inputs.

*   **Agent Integration:** The existing Coder, Auditor, and QA Analyst agents will be integrated as nodes in the LangGraph orchestration graph.

The workflow for this cycle will be managed by the LangGraph graph and will look as follows:
1. The graph is initiated with the cycle ID.
2. The Coder agent is invoked to write the initial version of the code.
3. The unit tests are run. If they fail, the graph transitions back to the Coder agent with the test results as feedback.
4. If the unit tests pass, the Auditor agents are invoked. If they find issues, the graph transitions back to the Coder agent with the audit report as feedback.
5. If the audit is clean, the QA Analyst agent is invoked. If the UATs fail, the graph transitions back to the Coder agent with the UAT results as feedback.
6. If the UATs pass, the graph transitions to a final "success" state, and the cycle is complete.

## 3. Design Architecture

The design architecture for Cycle 05 focuses on the implementation of the LangGraph orchestration engine.

**File Structure:**

The following new and updated files will be created within `dev_src/ac_cdd_core/`:

*   `orchestration/graph.py`: A new file to define the LangGraph orchestration graph.
*   `orchestration/nodes.py`: A new file to define the functions that will serve as the nodes in the graph (e.g., `run_coder`, `run_auditor`).
*   `services/coder.py`: The `CoderService` will be updated to use the new orchestration graph.

**Class and Function Definitions:**

*   **`orchestration/graph.py`:**
    *   `create_graph()`: A function that builds and returns the LangGraph `StatefulGraph` instance. It will define the graph's state schema and the edges between the nodes.

*   **`orchestration/nodes.py`:**
    *   `run_coder(state)`: A function that invokes the Coder agent.
    *   `run_unit_tests(state)`: A function that runs the unit tests in the sandbox.
    *   `run_auditor(state)`: A function that invokes the Auditor agent.
    *   `run_qa_analyst(state)`: A function that invokes the QA Analyst agent.

*   **`services/coder.py`:**
    *   `CoderService.run()`: This method will be updated to create an instance of the orchestration graph and then run it.

This design cleanly separates the orchestration logic from the agent and service implementations, making the system more modular and easier to understand.

## 4. Implementation Approach

The implementation of Cycle 05 will be carried out in the following steps:

1.  **Define the graph state:** The Pydantic model for the graph's state will be defined. This will include fields for the current code, test results, audit reports, and UAT results.

2.  **Implement the graph nodes:** The functions for each of the nodes in the graph will be implemented in `orchestration/nodes.py`.

3.  **Build the graph:** The `create_graph` function will be implemented in `orchestration/graph.py`. This will involve adding the nodes to the graph and defining the conditional edges that control the flow of execution.

4.  **Update `CoderService`:** The `CoderService` will be refactored to use the new orchestration graph.

5.  **Write unit tests:** Unit tests will be written for the orchestration graph and its nodes. The agents and services will be mocked to allow the graph's logic to be tested in isolation.

6.  **Write end-to-end tests:** End-to-end tests will be created for the `run-cycle` command. These tests will simulate a complete development cycle, including several iterations of the fixing loop, to verify that the orchestration logic is working correctly.

## 5. Test Strategy

The test strategy for Cycle 05 will focus on ensuring the correctness and robustness of the LangGraph-based orchestration engine.

**Unit Testing Approach:**

*   **Orchestration Graph:** The unit tests will focus on the logic of the graph itself. We will create a mock `StatefulGraph` and test the conditional edges to ensure that the graph transitions between states correctly based on the inputs. The functions in `orchestration/nodes.py` will be tested individually, with the agents and services they call being mocked.

**Integration and End-to-End Testing Approach:**

The end-to-end tests for the `run-cycle` command will be the primary method for validating the functionality of this cycle. These tests will be designed to cover various scenarios for the iterative fixing loop:

*   **Single-Pass Success:** A test where the code is correct on the first try and passes all the quality gates without needing any fixes.

*   **Fixing Unit Test Failures:** A test where the Coder agent's initial code fails the unit tests. The test will mock the agent to produce correct code on the second attempt and will verify that the graph correctly loops back to the coder and eventually succeeds.

*   **Fixing Audit Issues:** A test where the code passes the unit tests but fails the audit. The test will verify that the graph correctly identifies the audit failure and re-engages the Coder agent to fix the issues.

*   **Fixing UAT Failures:** A test where the code passes the unit tests and audit but fails the UATs. The test will verify that the graph correctly handles the UAT failure and loops back to the Coder agent.

*   **Maximum Iterations:** A test to verify that the graph correctly handles the case where the Coder agent is unable to fix the issues within a predefined maximum number of iterations, preventing an infinite loop.
