import math

import pytest
from scipy.integrate import quad
from scipy.stats import gamma

from fraser.model import (
    exponential_threshold,
    gamma_parameters,
    isolation_threshold,
)


@pytest.mark.parametrize("theta", [0.05, 0.3, 0.8, 0.95])
@pytest.mark.parametrize("variance_ratio", [0.1, 0.3, 0.5, 1.0])
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
    assert exponential_threshold(theta, isolation, tracing=0) == expected


@pytest.mark.parametrize("theta", [0.05, 0.3, 0.8])
def test_perfect_isolation_and_tracing(theta):
    assert exponential_threshold(theta, isolation=1, tracing=1) == pytest.approx(
        (2 - theta) / theta
    )


@pytest.mark.parametrize("theta", [0, 0.3, 1])
def test_no_interventions(theta):
    assert exponential_threshold(theta, isolation=0, tracing=0) == 1


@pytest.mark.parametrize("theta", [0.01, 0.1, 0.5, 0.99])
@pytest.mark.parametrize("isolation", [0.75, 0.9, 1])
def test_tracing_increases_analytical_threshold(theta, isolation):
    assert exponential_threshold(theta, isolation, tracing=0) >= isolation_threshold(
        theta, isolation
    )


def test_time_scale_invariance():
    shape, infectiousness_scale, incubation_scale = gamma_parameters(0.3, 0.1)
    actual = gamma_parameters(0.3, 0.1, generation_time=30)
    assert actual == pytest.approx(
        (shape, 10 * infectiousness_scale, 10 * incubation_scale)
    )


def test_endpoints():
    assert isolation_threshold(1, 0.75) == 1
    assert exponential_threshold(1, 0.9, 1) == 1
    assert isolation_threshold(0.2, 0) == 1
