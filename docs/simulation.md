# Simulation methods

## Model

This is a homogeneous branching process with no susceptible depletion.
Each run starts with many seed infections.
These infections are essentially independent simulations because contact tracing is forward-only.
The paper used 100 replicates to locate zero growth; here each candidate uses three many-seed runs.

Infectiousness is $\beta(\tau) = R_0g(\tau)$, with $g$ a gamma density of mean $T_g = 1$.
Incubation times are independent gamma draws with mean $T_s$ and the same shape $k$.
Their scales are $1/k$ and $T_s/k$, respectively.
The paper's low/high variance cases correspond to $k = 10$ and $k = 2$.
The presymptomatic transmission fraction follows from the two distributions:

$$
  \theta = I_{T_s/(1+T_s)}(k,k)
$$

where $I$ is the regularized incomplete beta function.
All infections eventually develop symptoms; all times use units of $T_g$.
To cover the same theta range for both shapes, choose a target theta grid and infer the symptom mean:

$$
  z = I^{-1}_{\theta}(k,k),\qquad T_s = \frac{z}{1 - z}.
$$

## Timestep update

1. Advance infection ages by $\delta$.
1. Process infections.
   People neither isolated nor quarantined generate Poisson infections with mean $\beta(\tau)\delta$.
   Record infector-infectee links.
   Draw symptom onset times.
1. Process symptoms.
   Newly symptomatic, non-quarantined people isolate with probability $\varepsilon_I$.
   Newly symptomatic, previously quarantined people are isolated with probability 1.
1. Process tracing.
   Newly isolated people quarantine their (non-isolated, non-quarantined) infectees independently with probability $\varepsilon_T$.
   If a newly quarantined person is already symptomatic, they in turn trace their infectees.
   Iterate.

Note that the paper's algorithm calls transmitting people "not isolated"; we interpret this to exclude quarantine, consistent with its mathematical model (page 6148) and intervention description (page 6150).

## Defaults

  | Setting                                  | Value                          |
  | ---------------------------------------- | ------------------------------ |
  | Initial infections per run               | 300                            |
  | Timestep / horizon / incidence bin width | 0.02 / 10 / 1                  |
  | Target theta values                      | 0.05, 0.2, 0.4, 0.6, 0.8, 0.95 |
  | Shapes                                   | 2, 10                          |
  | Isolation                                | 1, 0.9, 0.75                   |
  | Tracing efficacy                         | 1                              |
  | Population cap                           | $10^6$ cumulative infections   |

## Growth rate search

Bin new infections, excluding seeds, and fit log incidence against bin midpoint from the start.
Include the first zero as 1, then discard later bins; with no zero, use the full recorded history.

For each target theta, shape, and isolation efficacy, calibrate $T_s$ and search with tracing enabled:

1. Evaluate $R_0 = 1$ and $R_0 = 15$ using the median growth fit from three independent runs each.
2. Report out-of-range crossings or unavailable/inconsistent endpoint fits explicitly.
3. Bisect: positive median growth lowers the upper bound; otherwise raise the lower bound.
4. Stop at bracket width 0.5 and report its midpoint.
