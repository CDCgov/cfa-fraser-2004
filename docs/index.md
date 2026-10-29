# Re-implementation of Fraser et al. (2004)

See [Fraser et al. (2004)](https://doi.org/10.1073/pnas.0307506101) for motivation and explanation of the math.
This package implements the closed-form [deterministic equations](deterministic.md).
The [stochastic simulator](simulation.md) produces exploratory incidence timecourses with isolation and tracing.

## Results

![](figures/deterministic.svg)

The solid black lines reproduce Figure 3's isolation-only thresholds.
The dashed blue lines use corrected Equation 18, assuming independent isolation and tracing, exponential-distributed incubation, exponential-distributed infectiousness, and 100% tracing efficacy.
Gray rectangles show the paper's historical disease ranges.
