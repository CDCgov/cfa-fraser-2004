"""Deterministic isolation and independent-tracing thresholds."""

import math

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigvals
from scipy.special import betaincinv, gammaln, logsumexp, xlogy
from scipy.stats import gamma


def _probability(value: float, name: str) -> None:
    if not 0 <= value <= 1:
        raise ValueError(f"{name} must be between 0 and 1")


def isolation_threshold(theta: float, isolation: float = 1.0) -> float:
    """Return critical R0 from Equation 6; control requires R0 below it."""
    _probability(theta, "theta")
    _probability(isolation, "isolation")
    remaining = 1 - isolation + isolation * theta
    return 1 / remaining if remaining else math.inf


def exponential_threshold(
    theta: float, isolation: float = 1.0, tracing: float = 1.0
) -> float:
    """Return Equation 18's threshold with its denominator typo corrected."""
    _probability(theta, "theta")
    _probability(isolation, "isolation")
    _probability(tracing, "tracing")
    remaining = (
        (1 - isolation) * (1 - tracing)
        + (isolation * (1 - tracing) + (1 - isolation) * tracing) * theta
        + isolation * tracing * theta / (2 - theta)
    )
    return 1 / remaining if remaining else math.inf


def gamma_parameters(
    theta: float, variance_ratio: float, generation_time: float = 3.0
) -> tuple[int, float, float]:
    """Return common shape, infectiousness scale, and incubation scale.

    Both distributions have variance / mean**2 equal to variance_ratio.
    Integer gamma shapes are supported, including the paper's 10 and 2,
    and the exponential benchmark's 1. Time units are arbitrary.
    """
    if not 0 < theta < 1:
        raise ValueError("theta must be strictly between 0 and 1")
    if not math.isfinite(variance_ratio) or variance_ratio <= 0:
        raise ValueError("variance_ratio must be positive and finite")
    shape = round(1 / variance_ratio)
    if shape < 1 or not math.isclose(shape, 1 / variance_ratio):
        raise ValueError("variance_ratio must be the reciprocal of an integer")
    if not math.isfinite(generation_time) or generation_time <= 0:
        raise ValueError("generation_time must be positive and finite")
    fraction = betaincinv(shape, shape, theta)
    infectiousness_scale = generation_time / shape
    incubation_scale = infectiousness_scale * fraction / (1 - fraction)
    return shape, infectiousness_scale, incubation_scale


def _log_survival(age: np.ndarray, shape: int, scale: float) -> np.ndarray:
    scaled = np.asarray(age) / scale
    powers = np.arange(shape)
    return -scaled + logsumexp(
        xlogy(powers, scaled[..., None]) - gammaln(powers + 1), axis=-1
    )


def tracing_threshold(
    theta: float,
    isolation: float = 1.0,
    tracing: float = 1.0,
    variance_ratio: float = 0.1,
    generation_time: float = 3.0,
    quadrature_order: int = 96,
    tail_probability: float = 1e-10,
) -> float:
    """Return critical R0 from Equation 16, not the correlated simulation.

    Gauss-Legendre quadrature discretizes the next-generation operator.
    Its spectral radius at R0=1 determines the critical R0 by inversion.
    Isolation and tracing are independent in this approximation.
    """
    _probability(isolation, "isolation")
    _probability(tracing, "tracing")
    shape, infectiousness_scale, incubation_scale = gamma_parameters(
        theta, variance_ratio, generation_time
    )
    if tracing == 0:
        return isolation_threshold(theta, isolation)
    if not isinstance(quadrature_order, int) or quadrature_order < 2:
        raise ValueError("quadrature_order must be an integer of at least 2")
    if not 0 < tail_probability < 1:
        raise ValueError("tail_probability must be strictly between 0 and 1")
    upper = gamma.isf(tail_probability, shape, scale=infectiousness_scale)
    nodes, weights = leggauss(quadrature_order)
    ages = (nodes + 1) * upper / 2
    weights = weights * upper / 2
    log_survival = _log_survival(ages, shape, incubation_scale)
    conditional_survival = np.exp(
        _log_survival(ages[:, None] + ages[None, :], shape, incubation_scale)
        - log_survival[None, :]
    )
    infectiousness = gamma.pdf(ages, shape, scale=infectiousness_scale)
    isolation_factor = 1 - isolation + isolation * np.exp(log_survival)
    operator = (
        (infectiousness * isolation_factor)[:, None]
        * (1 - tracing + tracing * conditional_survival)
        * weights[None, :]
    )
    radius = float(np.max(np.abs(eigvals(operator))))
    return 1 / radius if radius else math.inf
