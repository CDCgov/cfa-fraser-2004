"""Deterministic isolation and independent-tracing thresholds."""

from scipy.special import betaincinv


def isolation_threshold(theta: float, isolation: float) -> float:
    """Critical R0 from Equation 6; control requires R0 below it."""
    return 1 / (1 - isolation + isolation * theta)


def exponential_threshold(theta: float, isolation: float, tracing: float) -> float:
    """Equation 18's threshold (with typo corrected)."""
    return 1 / (
        (1 - isolation) * (1 - tracing)
        + (isolation * (1 - tracing) + (1 - isolation) * tracing) * theta
        + isolation * tracing * theta / (2 - theta)
    )


def gamma_parameters(
    theta: float, variance_ratio: float, generation_time: float = 3.0
) -> tuple[float, float, float]:
    """Calibrate gamma parameters for stochastic simulations.

    Both distributions have variance / mean**2 equal to variance_ratio.
    Returns common shape, infectiousness scale, and incubation scale.
    Time units are arbitrary; no deterministic operator is evaluated.
    """

    shape = 1 / variance_ratio
    fraction = betaincinv(shape, shape, theta)
    infectiousness_scale = generation_time / shape
    incubation_scale = infectiousness_scale * fraction / (1 - fraction)
    return shape, infectiousness_scale, incubation_scale
