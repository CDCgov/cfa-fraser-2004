# Simulation protocol for Figure 3

**Specification for future implementation.** The current package computes deterministic thresholds only.
Proposed settings below are not historical values or validated defaults.

## Role of the simulation

Figure 3's black lines come from Equation 6.
Its colored lines were found by repeating the stochastic individual-based simulation 100 times for candidate $(R_0,\theta)$ combinations and finding where incidence neither grew nor declined on average.
Simulation also demonstrated the limitations of the independent tracing approximation in Equations 14-18.

In the realistic simulation, one symptom time determines both isolation and the opportunity to trace that person's contacts.
These events are correlated.
The paper also describes a modified simulation with independent clocks to diagnose the analytical approximation, not to construct the realistic curves.
It reports agreement at perfect isolation for the distributions tested and potentially substantial overestimation of tracing benefit at imperfect isolation.

## Model specified by the paper

  | Input                   | Choice                                                       |
  | ----------------------- | ------------------------------------------------------------ |
  | Population              | Early outbreak, homogeneous mixing, no susceptible depletion |
  | Infectiousness          | Common rate $\beta(\tau) = R_0g(\tau)$                       |
  | Incubation              | Independent individual draws with survival function $S$      |
  | Distributions           | Gamma infectiousness and incubation densities                |
  | Variance / mean squared | 0.1 or 0.5 for both distributions                            |
  | Isolation efficacy      | 1, 0.9, or 0.75                                              |
  | Tracing efficacy        | 1 for colored curves, 0 for isolation-only checks            |
  | Explicit delays         | None                                                         |
  | Replicates              | 100 per evaluated parameter combination                      |
  | Step size               | Reduced until convergence; final value not given             |

Use the [same gamma calibration](deterministic.md#gamma-calibration) as in the deterministic calculation.
This is a branching process, not a finite-population SIR model.
Do not remove infectious people at an arbitrary recovery age without controlling the omitted gamma-tail mass.

## State and timestep algorithm

Store each person's infection time, symptom time, infector, infected children, isolation-success indicator, isolation/quarantine status, and traced-once flag.
Retain removed individuals: future symptom events can still initiate tracing.

At time zero, infect $Y_i$ individuals at age zero and draw incubation durations and independent Bernoulli isolation-success indicators with probability $\varepsilon_I$.
Use an explicit NumPy random generator.

Follow the paper's stated timestep ordering:

1. Advance infection ages by $\delta$.
2. Each person neither isolated nor quarantined generates $\mathrm{Poisson}(\beta(\tau)\delta)$ infections at the updated age.
3. Add newborn infections at age zero; draw their incubation durations and isolation indicators, and record parent-child links.
   Newborns do not transmit during this same step.
4. Process symptoms and removals.
   A newly symptomatic person is isolated if isolation succeeds.
   An isolated symptomatic person has each infected child quarantined independently with probability $\varepsilon_T$.

Quarantine stops transmission immediately.
A quarantined presymptomatic contact does **not** trigger tracing until symptoms appear.
Only forward tracing along infected parent-child links is included: no backward tracing, uninfected-contact screening, or recursive presymptomatic tracing.
Isolation is all-or-nothing, not a partial rate reduction for every symptomatic person.
These choices give the same isolation-only mean but need not be interchangeable once tracing is introduced.

### Event conventions not fully specified

Proposed choices that must be recorded and checked:

- A quarantined person is already removed.
  At their symptoms they can trace children even if their independent isolation indicator was false.
  This interprets a quarantined symptomatic person as an isolated case.
  Check the alternative interpretation for sensitivity at imperfect isolation.
- If a quarantined child is already symptomatic, process their tracing in the same removal phase.
  Use a queue and traced-once flags.
  Cascades may pass through symptomatic people, never through presymptomatic people.
- Symptoms occurring within a step are processed after its transmissions, as listed in the paper.
  This artificial delay is at most $\delta$ and must vanish under timestep refinement.
- Interventions operate from time zero.
  Seeds have no recorded infector and cannot be traced from an external source, but can initiate tracing themselves.
  Exclude the initialization transient when estimating growth.

These are reproducible interpretations, not claims about the authors' code.
If they change the limiting curves materially, report the ambiguity rather than silently tuning conventions to match the image.

## Omitted settings and proposed starting values

Time is expressed in generation-time units; $T_g = 3$ arbitrary units fixes the overall scale without changing the idealized threshold.

  | Setting omitted or incompletely specified | Proposed starting choice                                                                       |
  | ----------------------------------------- | ---------------------------------------------------------------------------------------------- |
  | Initial infections $Y_i$                  | 1,000 per replicate                                                                            |
  | Numerical timestep                        | $T_g/100$; repeat with $T_g/200$ and $T_g/400$                                                 |
  | Horizon                                   | $40T_g$; lengthen to test stability                                                            |
  | Growth fitting window                     | $20T_g$ to $40T_g$                                                                             |
  | Incidence bins                            | Width $T_g$, counting new infections                                                           |
  | Theta grid                                | 0.02 to 0.98 by 0.02, refining curved regions                                                  |
  | R0 bracket                                | Start at 1 and 2, doubling the upper bound until growth is established                         |
  | R0 search safeguard                       | 128; if unbracketed, report unresolved or above bound, not threshold 128                       |
  | Precision target                          | R0 bracket width 0.02, subject to Monte Carlo uncertainty                                      |
  | Replicates                                | Initially 100; increase near ambiguous boundaries                                              |
  | Randomness                                | Master seed 2004 and recorded independent `SeedSequence` streams per evaluation/replicate      |
  | Resource safeguard                        | Initially $10^6$ infections per replicate; reaching it invalidates a complete-history estimate |
  | Uncertainty                               | Bootstrap entire independent replicate histories; initially 1,000 resamples                    |

Resource limits are not epidemic outcomes.
Do not classify capped runs as extinct, drop them from the mean, or fit truncated histories as complete.
Increase resources or report the evaluation unresolved.
The plot's maximum R0 of 15 must not constrain the numerical search.

## Growth estimator and threshold search

The paper does not report its exact estimator.
A proposed reproducible method:

1. Simulate independent replicates at each parameter combination.
2. Compute binned **mean incidence including extinct runs' zeros**.
   Do not condition on survival.
3. Fit log mean incidence against bin midpoint over the late-time window.
   Its slope estimates exponential growth rate $r$.
   Do not add pseudocounts.
   If a required mean bin is zero, the slope is unresolved; increase sampling or reconsider the window.
   All-replicate extinction alone does not establish subcriticality.
4. Bootstrap whole replicate histories and repeat the fit.
   Record failures caused by zero bins rather than silently discarding bootstrap samples.
5. Bracket negative and positive growth.
   Bisect only when the sign is adequately resolved.
   If uncertainty spans zero, increase sampling rather than trusting a noisy sign.
6. Report the R0 bracket and Monte Carlo uncertainty.
   Repeat with longer horizons and different fitting windows to check asymptotic behavior.

The target is $r = 0$, not a chosen extinction probability.
Even supercritical finite outbreaks can die out; near criticality a few large outbreaks may dominate the mean.
More seeds, replicates, and horizon checks address distinct problems and are not substitutes for one another.

## Required validation

1. With tracing off, recover Equation 6 for both variances and all three isolation efficacies within Monte Carlo uncertainty.
2. With no interventions, recover a mean of $R_0$ lifetime infections per case, allowing for numerical error in the timestep sum of $\beta(\tau)\delta$.
3. With zero isolation efficacy and no initially quarantined people, tracing must have no effect: nobody can initiate it.
   The independent analytical approximation need not satisfy this realistic-model check.
4. With perfect isolation, compare the realistic simulation's thresholds to Equation 16 for shapes 1, 2, and 10, as checked by the paper.
5. In a separate diagnostic mode, draw independent isolation and tracing clocks, and let the tracing clock act independently of successful isolation.
   Compare with the independent analytical model.
   Do not label diagnostic curves as realistic Figure 3 results.
6. Refine timestep, lengthen horizon, vary fitting windows, and increase seeds and replicates.
   Require changes small compared with reported uncertainty; check ambiguous same-step conventions separately.
7. Test a hand-built tree: quarantining a presymptomatic child blocks its transmission but does not trace its children until symptoms.
   Each eligible person traces at most once.

Only after these checks should simulated curves replace or accompany the current deterministic upper bounds.
