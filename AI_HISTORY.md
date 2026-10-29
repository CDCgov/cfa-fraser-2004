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

### Singular endpoint simplification

- **Task completed**: Removed the test requiring a division-by-zero guard to return infinity, and removed the deterministic plot's special theta-zero evaluation.
  Retained the simple formulas and substantive SVG tests.

### Rough growth-rate grid

- **Task completed**: Added a simple log-incidence slope fit from the first bin through the first zero, replacing that zero with 1 and ignoring later bins.
  Expanded the editable experiment to a 4-by-4 R0/symptom-mean grid for two gamma shapes and two policies.
  Generated an annotated Altair growth-rate SVG with population-limit markers.
  Kept seed infections separate from plotted incidence, documented the rough fitting rule, and added fitting/grid tests.

### Critical-R0 search

- **Task completed**: Added bisection over R0=1-15 at fixed symptom-onset mean and gamma shape.
  Each candidate uses the median rough growth fit from three independent seeds; searches stop at a bracket width of 0.25.
  Added explicit out-of-range and unresolved statuses, the `critical` command, editable search settings, and an Altair curve with numerical brackets (not confidence intervals).
  Ran all 16 default searches and added search, command, and SVG tests.

### Single-figure workflow

- **Task completed**: Removed timecourse, growth-grid, and standalone deterministic plots and their command paths.
  The default command now generates one simulation-led control figure with faint dashed analytical overlays and historical disease tiles.
  Replaced the simulation walkthrough with concise methods, removed code examples, and updated documentation and plot tests.

### Efficacy panels and theta coverage

- **Task completed**: Switched the experiment to target theta values from 0.05 to 0.95, with inverse gamma calibration for symptom means and the forward computation retained in the simulator.
  Panels now show isolation efficacies 100%, 90%, and 75%; blue/red distinguish simulated shapes 10/2.
  Kept Equation 6 as a black analytical curve and Equation 18 as a dashed gray reference.
  Eliminated isolation-only simulations and adopted faster preview settings: 300 seeds, timestep 0.02, and bracket tolerance 0.5.
  Updated methods and tests.
