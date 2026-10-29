# AI History

## 2026-10-29

- **Task completed**: Implemented deterministic outbreak-control thresholds from Fraser et al. (2004), including the isolation-only curves and Equation 16 tracing approximation.
  Generated minimal Altair SVGs with disease-range rectangles, dashed tracing bounds, and an R0 axis from 0 to 15 with ticks every 3.
  Documented the core mathematics, Equation 18 typo, and future simulation protocol.
  Validated with 65 passing tests, lint, and a strict documentation build.
  Created this history file.
- **Model**: GPT-6.1 Sol
- **Agent environment**: GitHub Copilot Chat in VS Code on Linux.

### Closed-form simplification

- **Task completed**: Removed the Equation 16 quadrature, eigensolver, and numerical benchmark plot.
  The deterministic SVG now uses Equation 6 and corrected Equation 18, with one dashed exponential tracing curve per panel.
  Replaced the numerical Equation 18 check with algebraic special-case tests.
  Documentation starts from Equation 17, derives the closed form, and retains the denominator typo note.
  Gamma calibration and sampling are reserved for the future stochastic simulation; disease rectangles and the 0-15 axis with ticks every 3 are preserved.

### Timestep simulation

- **Task completed**: Implemented one-population, many-seed branching simulations with vectorized Poisson transmission, gamma symptom draws, and symptom-triggered forward tracing.
  Added a parameter-grid command producing linear/log Altair incidence SVGs and explicit population-limit reporting.
  Documented event ordering, quarantine conventions, defaults, and the distinction between seed population and incidence.
  Recorded full-window growth fitting as a deferred next step.
  Added tracing-tree, renewal-mean, reproducibility, stopping, command, and SVG tests.
  Also restored the infinite isolation-threshold endpoint needed by the existing plot.

### CLI simplification

- **Task completed**: Moved the simulation grid and numerical parameters into a single editable `SIMULATION_SETTINGS` block in the package entry point.
  Removed simulation parameter flags, keeping the `simulate` command and `--output-dir` option.
  Updated the command tests and documentation; the Python simulation API is unchanged.
