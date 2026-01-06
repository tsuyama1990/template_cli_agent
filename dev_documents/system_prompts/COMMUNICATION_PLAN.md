# Communication Plan: MLIP-AutoPipe

This document outlines the communication strategy for the MLIP-AutoPipe project to ensure all stakeholders are informed and engaged throughout the development lifecycle.

## 1. Stakeholders

*   **Project Lead/Owner:** The primary decision-maker for the project.
*   **Development Team:** The engineers responsible for designing, building, and testing the software.
*   **Scientific Advisors/End-Users:** Domain experts and target users who provide scientific requirements and validate the utility of the generated data.

## 2. Communication Channels & Cadence

| Channel | Purpose | Frequency | Audience |
|---|---|---|---|
| **Daily Stand-up** | Quick synchronization on progress, daily goals, and blockers. | Daily | Development Team |
| **Cycle Planning Meeting** | Review backlog, define sprint goals, and plan tasks for the upcoming cycle. **Includes a standing agenda item for Risk Review.** | Start of each 2-week cycle | Development Team, Project Lead |
| **End-of-Cycle Demo & Review** | Demonstrate the features completed in the cycle, gather feedback, and review against the cycle goals. | End of each 2-week cycle | All Stakeholders |
| **GitHub** | Source code management, issue tracking, and asynchronous technical discussions via Pull Requests. | Continuous | Development Team |
| **Shared Project Board** | Kanban-style board (e.g., GitHub Projects) to track the status of all tasks. | Continuous | All Stakeholders |

## 3. Reporting

*   **Cycle Plan:** At the beginning of each cycle, the plan and goals for the sprint will be shared with all stakeholders.
*   **Weekly Summary:** A brief weekly email summary from the Project Lead highlighting key accomplishments, progress, and any significant risks.
*   **Cycle Retrospective:** A private meeting for the Development Team after the End-of-Cycle Demo to discuss internal process improvements.

## 4. Stakeholder Feedback Incorporation (Operational Plan)

A structured feedback loop is critical to ensuring the project's output is scientifically valid and useful.

**Process:**

1.  **Demonstration with Scientific Context:** At the End-of-Cycle Demo, the Development Team will present the new functionality. **This presentation will be centered around the UAT Jupyter Notebook for that cycle**, which allows for a live demonstration of the tool on a real scientific problem.
2.  **Targeted Feedback Collection:** During the demo, feedback will be actively solicited on:
    *   The **scientific validity** of the generated structures.
    *   The **diversity** of the dataset.
    *   The usability of the CLI and configuration files.
    *   Any discrepancies between the output and the expected physical behavior.
3.  **Triage and Prioritization:** The Project Lead and Development Team will hold a follow-up meeting to triage the collected feedback. Each item will be categorized (e.g., Scientific Bug, Usability Issue, New Feature Request) and prioritized based on its impact on the project's scientific goals.
4.  **Backlog Integration:** High-priority items will be converted into well-defined issues in the GitHub backlog. These issues will be considered for inclusion in the next cycle's planning meeting.
5.  **Transparent Communication:** The Project Lead will share the outcome of the triage meeting with stakeholders, explaining which feedback items have been incorporated into the next cycle's plan and providing a rationale for any items that are deferred. This ensures a transparent and collaborative process.
