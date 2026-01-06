# Risk Management Plan: MLIP-AutoPipe

This document outlines the potential risks for the MLIP-AutoPipe project and the strategies to mitigate them. **This document will be reviewed and updated at the beginning of each development cycle planning session.**

## 1. Technical Risks

| Risk ID | Description | Probability | Impact | Mitigation Strategy | Contingency Plan |
|---|---|---|---|---|---|
| T-01 | **MLIP Model Integration Issues:** The MACE/SevenNet models are complex external dependencies. Updates could introduce breaking API changes, or bugs could affect simulation stability. | Medium | High | - Lock library versions in `pyproject.toml` using `uv.lock`.<br>- Develop a dedicated integration test for the MLIP calculator to validate its basic functionality.<br>- Abstract the calculator interaction behind a dedicated wrapper class. | Revert to the last known working version and open an issue with the model's developers, providing a minimal reproducible example. |
| T-02 | **Performance Bottlenecks:** The parallel MD simulation may not scale as expected due to the GIL, data serialization overhead, or I/O contention. | Medium | Medium | - Use `ProcessPoolExecutor` for CPU-bound tasks.<br>- Implement "late-binding" of the calculator within each worker process to avoid pickling large models.<br>- Ensure each worker process writes to a unique temporary file. | Profile the application to identify the exact bottleneck. Explore alternative parallelization libraries or more advanced inter-process communication if needed. |
| T-03 | **Physics Violation Explosion:** High-temperature simulations are chaotic. A significant portion of simulations might fail due to atoms getting too close, leading to wasted computation. | High | Medium | - **Integrate the ZBL potential for short-range repulsion as a core feature.**<br>- Implement robust error handling within each simulation worker to catch physics violations, log the failing structure for analysis, and exit gracefully. | Implement a "dynamic temperature ramp," starting simulations at a lower temperature and gradually increasing it to allow the system to relax gently. |

## 2. Schedule and Project Management Risks

| Risk ID | Description | Probability | Impact | Mitigation Strategy | Contingency Plan |
|---|---|---|---|---|---|
| S-01 | **Underestimation of Scientific Complexity:** The complexity of implementing advanced scientific features (e.g., Hybrid MD/MC, FPS sampling) may be greater than anticipated, leading to sprint delays. | Medium | High | - Break down large scientific tasks into smaller, verifiable sub-tasks.<br>- Prioritize the core functionality and define a clear MVP for each cycle.<br>- Foster a culture of early and open communication about implementation challenges. | Re-evaluate the feature's priority in consultation with scientific advisors. We may decide to move it to a subsequent cycle or simplify its initial implementation (e.g., use a simpler diversity metric than FPS). |
| S-02 | **Dependency on Key Personnel:** In a small team, the unavailability of a key person could significantly impact the timeline. | Low | High | - Enforce a strict code review policy to ensure knowledge sharing.<br>- Maintain high-quality, up-to-date documentation for all components.<br>- Use a shared project board for full visibility of task status. | Conduct a rapid knowledge transfer session. Adjust the project scope or timeline in consultation with stakeholders. |

## 3. External Risks

| Risk ID | Description | Probability | Impact | Mitigation Strategy | Contingency Plan |
|---|---|---|---|---|---|
| E-01 | **Upstream Library Bugs:** A critical bug in a core dependency (e.g., ASE, `mace-torch`) could block development or produce scientifically incorrect results. | Low | High | - Pin dependency versions in `uv.lock`.<br>- Our comprehensive integration test suite will help catch regressions introduced by dependencies. | First, attempt to find a workaround. If none is available, report the bug to the library's maintainers with a minimal reproducible example. We can also investigate pinning to an older, stable version. |
