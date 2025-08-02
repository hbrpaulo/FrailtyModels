# AGENTS.md

## Workflow
- Keep pull requests (PRs) small and focused on one task group.
- Merge each PR before starting the next one to avoid conflicts.
- Use feature branches named after the work item (e.g., `feature/split-funcoes`).

---

## Code Organization
- Legacy R code should be modularized into `scripts - R/legacy/{analysis,interpretation,report}`.
- Each R function lives in its own file, uses `snake_case`, and includes
`roxygen2` documentation.
- Place example notebooks or R Markdown files under `usage_examples/legacy`.

---

## Coding Style
- **R:** `roxygen2` headers, one function per file, descriptive variable names.
- **Python:** Follow PEP 8 style; keep lines ≤79 characters and write
docstrings for public functions.

---

## Testing
- **Python tests:** run `pytest` at repository root.
- **R tests:** run `R -q -e "testthat::test_dir('tests/legacy')"` after modifying any R code.

---

## Commit & PR Guidelines
- Commit messages should be in the imperative mood (e.g., “Add density plot helper”).
- PR descriptions must include:
- **Summary:** brief description of the change.
- **Testing:** commands executed and their results.
- Reference any related issues or task IDs in the PR body.

---

## Documentation
- Update `README.md` and relevant `usage_examples/` when introducing new features.
- Include comments or docstrings explaining non-obvious code sections.