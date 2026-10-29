import math

import pytest
from scipy.integrate import quad
from scipy.stats import gamma

from fraser.model import (
    exponential_threshold,
    gamma_parameters,
    isolation_threshold,
    tracing_threshold,
)


@pytest.mark.parametrize("theta", [0.05, 0.3, 0.8, 0.95])
@pytest.mark.parametrize("variance_ratio", [0.1, 0.5, 1.0])
def test_calibrated_theta(theta, variance_ratio):
    shape, infectiousness_scale, incubation_scale = gamma_parameters(
        theta, variance_ratio
    )
    actual, _ = quad(
        lambda age: (
            gamma.pdf(age, shape, scale=infectiousness_scale)
            * gamma.sf(age, shape, scale=incubation_scale)
        ),
        0,
        math.inf,
    )
    assert actual == pytest.approx(theta, abs=1e-9)


@pytest.mark.parametrize("isolation", [0, 0.75, 0.9, 1])
@pytest.mark.parametrize("theta", [0.05, 0.4, 0.95])
def test_no_tracing(theta, isolation):
    expected = isolation_threshold(theta, isolation)
    assert tracing_threshold(theta, isolation, tracing=0) == expected
    assert exponential_threshold(theta, isolation, tracing=0) == expected


@pytest.mark.parametrize("theta", [0.05, 0.3, 0.8])
@pytest.mark.parametrize("isolation,tracing", [(1, 1), (0.9, 1), (0.75, 0.6)])
def test_exponential_operator(theta, isolation, tracing):
    actual = tracing_threshold(theta, isolation, tracing, variance_ratio=1)
    expected = exponential_threshold(theta, isolation, tracing)
    assert actual == pytest.approx(expected, rel=1e-7)


@pytest.mark.parametrize("variance_ratio", [0.1, 0.5])
@pytest.mark.parametrize("theta", [0.01, 0.1, 0.5, 0.99])
@pytest.mark.parametrize("isolation", [0.75, 0.9, 1])
def test_quadrature_convergence(theta, isolation, variance_ratio):
    coarse = tracing_threshold(theta, isolation, variance_ratio=variance_ratio)
    fine = tracing_threshold(
        theta, isolation, variance_ratio=variance_ratio, quadrature_order=160
    )
    assert coarse == pytest.approx(fine, rel=1e-4)
    assert coarse >= isolation_threshold(theta, isolation) * (1 - 1e-8)


def test_time_scale_invariance():
    expected = tracing_threshold(0.3, 0.9)
    assert tracing_threshold(0.3, 0.9, generation_time=30) == pytest.approx(expected)


def test_endpoints():
    assert isolation_threshold(0, 1) == math.inf
    assert isolation_threshold(1, 0.75) == 1
    assert exponential_threshold(1, 0.9, 1) == 1
    assert isolation_threshold(0.2, 0) == 1


@pytest.mark.parametrize("theta", [-0.1, 1.1, math.nan])
def test_invalid_probability(theta):
    with pytest.raises(ValueError):
        isolation_threshold(theta)


def test_invalid_gamma_shape():
    with pytest.raises(ValueError, match="reciprocal"):
        gamma_parameters(0.3, 0.3)
