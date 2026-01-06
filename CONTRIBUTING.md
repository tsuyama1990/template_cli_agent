# Contributing to MLIP-AutoPipe

We welcome contributions to the MLIP-AutoPipe project. To ensure a high standard of code quality, scientific validity, and maintainability, we ask that all contributors adhere to the following guidelines.

## 1. Development Workflow

All development work, from bug fixes to new features, should follow this structured workflow:

1.  **Create an Issue:** Before starting work, create a GitHub Issue that describes the bug or feature.
2.  **Fork and Branch:** Fork the repository and create a new feature branch from `main`. Branch names should be descriptive, e.g., `feature/add-ionic-generator`.
3.  **Develop and Test:** Write your code, ensuring you add or update tests that cover your changes.
4.  **Adhere to Quality Standards:** Before committing, ensure your code meets our quality standards by running the pre-commit hooks.
5.  **Submit a Pull Request (PR):** Submit a Pull Request to the main repository, referencing the original GitHub Issue.

## 2. Code Quality and Pre-commit Hooks

We enforce a strict set of code quality standards using automated tools. All contributors must use the provided pre-commit hooks.

**Setup:**

```bash
# Install pre-commit if you haven't already
pip install pre-commit
# Set up the git hooks
pre-commit install
```

On every `git commit`, the following checks will be run: `ruff format`, `ruff check`, and `mypy`. **A commit will be rejected if any of these checks fail.**

## 3. Code Review Process

All Pull Requests (PRs) must be reviewed and approved by at least one other member of the development team before they can be merged.

### **Standard Review Criteria:**

1.  **Correctness:** Does the code achieve the goals outlined in the associated issue?
2.  **Design Adherence:** Does the code align with the architecture in `SYSTEM_ARCHITECTURE.md`?
3.  **Test Coverage:** Is the code adequately covered by unit and integration tests?
4.  **Readability and Maintainability:** Is the code clean, well-organized, and easy to understand?
5.  **Documentation:** Is the code documented with docstrings? Have relevant project documents been updated?
6.  **CI Checks:** All automated checks must be passing.

### **Scientific Code Review Criteria:**

In addition to the standard criteria, all contributions that touch the scientific logic of the pipeline **must** be reviewed against these specific requirements:

7.  **Adherence to Physical Constraints:**
    *   Does the code correctly enforce physical rules? (e.g., The `AlloyGenerator` must prevent atomic overlap; the `IonicGenerator` must ensure charge neutrality).
    *   Reviewers should critically assess whether any changes could lead to physically implausible structures.

8.  **Correctness of Scientific Algorithms:**
    *   Is the implementation of scientific methods correct? (e.g., The logic for Monte Carlo swaps, the distance calculations in Farthest Point Sampling, the mixing of ZBL and MLIP potentials).
    *   Code in the `explorers/` and `sampling/` directories requires particularly rigorous scrutiny.

9.  **Performance Considerations:**
    *   Does the code avoid obvious performance pitfalls, especially in computationally intensive areas like the `md_engine`?
    *   Changes that could impact the performance of MD simulations must be carefully evaluated.

## 4. Definition of Done

A task is considered "Done" only when it meets all of the following criteria:

*   The code has been successfully implemented.
*   The code is covered by a comprehensive suite of passing unit and integration tests.
*   The code adheres to all quality standards and passes all pre-commit checks.
*   All relevant documentation (docstrings, `SPEC.md`, etc.) has been updated.
*   The changes have been approved in a code review against **both standard and scientific criteria**.
*   The Pull Request has been merged into the `main` branch.
