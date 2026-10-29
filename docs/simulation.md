# Simulation methods

Run `uv run fraser` to generate the combined [control figure](figures/control.svg).
Edit `SIMULATION_SETTINGS` and `SEARCH_SETTINGS` in the package entry point to change the experiment.
Only the output directory is a CLI option.

## Model

This is a homogeneous branching process with no susceptible depletion.
Each run starts with many independent seed infections; tracing couples people within a tree, not different trees.
The paper used 100 replicates to locate zero growth; here each candidate uses three many-seed runs.

Infectiousness is $\beta(\tau) = R_0g(\tau)$, with $g$ a gamma density of mean $T_g = 1$.
Incubation times are independent gamma draws with mean $T_s$ and the same shape $k$.
Their scales are $1/k$ and $T_s/k$, respectively.
The paper's low/high variance cases correspond to $k = 10$ and $k = 2$.
The presymptomatic transmission fraction follows from the two distributions:

$$
  heta = P(G < T) = I_{T_s/(1+T_s)}(k,k).
$$

Here $I$ is the regularized incomplete beta function.
All infections eventually develop symptoms; all times use units of $T_g$.
To cover the same theta range for both shapes, choose a target theta grid and infer the symptom mean:

$$
  z = I^{-1}_{\theta}(k,k),\qquad T_s = \frac{z}{1 - z}.
$$

`gamma_parameters` performs the inverse calibration with $T_g = 1$.
The simulator still accepts $T_s$ and computes theta in the forward direction with `presymptomatic_fraction`.

## Timestep update

1. Advance infection ages by $\delta$.
2. People neither isolated nor quarantined generate Poisson infections with mean $\beta(\tau)\delta$.
3. Add newborns, record parent-child links, and draw symptom times, rounded up to a timestep.
   Newborns first transmit on the next step.
4. At symptoms, isolate with probability $\varepsilon_I$; previously quarantined people are isolated with certainty.
   Isolated symptomatic cases quarantine eligible infectees independently with probability $\varepsilon_T$.
5. Newly quarantined, already-symptomatic contacts can immediately trace their children.
   Presymptomatic contacts wait until symptoms; each case traces at most once.

Quarantine blocks transmission but preserves later symptoms and tracing.
The paper's algorithm calls transmitting people "not isolated"; we interpret this to exclude quarantine, consistent with its mathematical model (p. 6148) and intervention description (p.
6150).
Certain isolation of quarantined symptomatic people and same-step symptomatic cascades are explicit implementation conventions.
Failed isolation is not retried; already removed contacts are unaffected by repeat tracing.
There is no backward tracing or presymptomatic recursive tracing.
Transmission precedes removal, so timestep refinement remains important.

## Defaults

  | Setting                                  | Value                                         |
  | ---------------------------------------- | --------------------------------------------- |
  | Initial infections per run               | 300                                           |
  | Timestep / horizon / incidence bin width | 0.02 / 10 / 1                                 |
  | Target theta values                      | 0.05, 0.2, 0.4, 0.6, 0.8, 0.95                |
  | Shapes                                   | 2, 10                                         |
  | Isolation / tracing efficacy             | 1, 0.9, 0.75 / tracing 1                      |
  | Population cap                           | $10^6$ cumulative infections, including seeds |
  | Master seed                              | 2004, with independent reproducible streams   |

The cap rejects the entire step that would exceed it; only completed steps are returned.
This is resource truncation, not extinction.
Zero-incidence timesteps do not stop the simulation.

## Growth and search

Bin new infections, excluding seeds, and fit log incidence against bin midpoint from the start.
Include the first zero as 1, then discard later bins; with no zero, use the full recorded history.
A final partial bin uses its actual duration.
Fewer than two bins makes the fit unavailable.
Finite fits from capped histories are retained.

For each target theta, shape, and isolation efficacy, calibrate $T_s$ and search with tracing enabled:

1. Evaluate R0=1 and 15 using the median growth fit from three independent runs each.
2. Report out-of-range crossings or unavailable/inconsistent endpoint fits explicitly.
3. Bisect: positive median growth lowers the upper bound; otherwise raise the lower bound.
4. Stop at bracket width 0.5 and report its midpoint.

Five midpoint evaluations give width 0.4375.
Tested R0 values and individual fits are retained in `CriticalResult.evaluations`.
This simple search reduces outlier sensitivity but can still choose the wrong half under noise.
Isolation-only thresholds use Equation 6 directly; they are not simulated.
Fewer seeds, a larger timestep, and a coarser bracket make these preview settings faster, at the cost of precision.

## Interpretation

Panels show isolation efficacy 100%, 90%, and 75%.
Black is the exact isolation-only threshold from Equation 6.
Blue and red curves with points show simulated tracing thresholds for shapes 10 and 2, respectively.
Vertical bars are **search brackets, not confidence intervals**.
Triangles mark crossings outside the search range; unresolved searches are omitted and counted.
The dashed gray curve is corrected [Equation 18](deterministic.md), evaluated at each panel's isolation efficacy and 100% tracing.
Equation 6 applies to either gamma shape; Equation 18 assumes exponential distributions and is a reference, not a gamma-model bound.
Gray rectangles show historical disease ranges, not probability distributions.

These are rough finite-horizon estimates, not a completed reproduction of Figure 3.
Check sensitivity to timestep, horizon, seed count, and random seed before interpreting small differences.
Near theta=1, symptoms become late and the fixed horizon may need extending.
