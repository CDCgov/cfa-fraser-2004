import numpy as np
import pytest
from scipy.stats import gamma

from fraser.simulation import (
    _People,
    _process_removals,
    presymptomatic_fraction,
    simulate,
)


@pytest.mark.parametrize("shape", [1, 2, 10])
def test_theta_from_means(shape):
    assert presymptomatic_fraction(1, shape) == pytest.approx(0.5)
    assert presymptomatic_fraction(0.5, shape) < 0.5
    assert presymptomatic_fraction(2, shape) > 0.5


def test_presymptomatic_quarantine_waits_for_symptoms():
    people = _People(1)
    people.add(np.array([-1]), 0, np.array([1]))
    people.add(np.array([0]), 0, np.array([3]))
    people.add(np.array([1]), 0, np.array([4]))
    rng = np.random.default_rng(1)
    assert _process_removals(people, 1, 1, 1, rng) == (1, 1)
    assert people.quarantined[1]
    assert not people.contacts_traced[1]
    assert not people.quarantined[2]
    assert _process_removals(people, 3, 0, 1, rng) == (1, 1)
    assert people.contacts_traced[1]
    assert people.quarantined[2]
    assert _process_removals(people, 3, 1, 1, rng) == (0, 0)


def test_symptomatic_quarantine_cascades():
    people = _People(4)
    people.add(np.array([-1]), 0, np.array([3]))
    people.add(np.array([0]), 0, np.array([1]))
    people.add(np.array([1]), 0, np.array([5]))
    rng = np.random.default_rng(2)
    assert _process_removals(people, 1, 0, 1, rng) == (0, 0)
    assert _process_removals(people, 3, 1, 1, rng) == (2, 2)
    assert people.quarantined[2]


def test_child_links_survive_growth_and_multiple_birth_batches():
    people = _People(1)
    people.add(np.array([-1, -1]), 0, np.array([3, 3]))
    people.add(np.array([0, 0, 1]), 1, np.array([5, 5, 5]))
    people.add(np.array([0, 1, 1]), 2, np.array([5, 5, 5]))
    assert _process_removals(people, 3, 1, 1, np.random.default_rng(1)) == (2, 6)
    assert people.quarantined[2:8].all()


def test_no_isolation_means_tracing_cannot_start():
    kwargs = dict(
        r0=1.2,
        symptom_mean=0.7,
        shape=2,
        isolation=0,
        initial_infections=100,
        duration=3,
        seed=17,
    )
    without = simulate(**kwargs, tracing=0)
    with_tracing = simulate(**kwargs, tracing=1)
    np.testing.assert_array_equal(without.incidence, with_tracing.incidence)
    assert not with_tracing.quarantines.any()


def test_first_step_poisson_mean_and_no_same_step_newborn_transmission():
    initial = 100_000
    result = simulate(
        2,
        1,
        2,
        isolation=0,
        initial_infections=initial,
        timestep=0.1,
        duration=0.1,
        seed=9,
    )
    expected = initial * 2 * gamma.pdf(0.1, 2, scale=0.5) * 0.1
    assert abs(result.incidence[1] - expected) < 6 * expected**0.5
    assert result.total_infections == initial + result.incidence[1]


def test_zero_incidence_step_does_not_stop_simulation():
    result = simulate(0, 0.1, 10, duration=1)
    assert result.completed
    assert result.times[-1] == 1
    assert not result.incidence.any()
    assert result.isolations.sum() == result.initial_infections


def test_population_limit_reports_only_complete_steps():
    result = simulate(
        20, 1, 2, initial_infections=100, timestep=0.1, duration=1, max_infections=100
    )
    assert not result.completed
    assert result.times.tolist() == [0]
    assert result.total_infections == 100


def test_reproducibility_and_accounting():
    first = simulate(1.5, 1, 2, initial_infections=100, duration=2, seed=20)
    second = simulate(1.5, 1, 2, initial_infections=100, duration=2, seed=20)
    np.testing.assert_array_equal(first.incidence, second.incidence)
    assert first.total_infections == first.initial_infections + first.incidence.sum()
    assert len(first.times) == len(first.incidence) == len(first.quarantines)


@pytest.mark.parametrize("shape", [2, 10])
@pytest.mark.parametrize("isolation", [0.75, 1])
def test_isolation_only_mean_matches_discrete_renewal(shape, isolation):
    initial = 10_000
    timestep = 0.02
    result = simulate(
        1.5,
        1,
        shape,
        isolation=isolation,
        tracing=0,
        initial_infections=initial,
        timestep=timestep,
        duration=2,
        seed=13,
    )
    ages = result.times[1:]
    kernel = (
        1.5
        * gamma.pdf(ages, shape, scale=1 / shape)
        * timestep
        * (
            1
            - isolation
            + isolation * gamma.sf(ages - timestep, shape, scale=1 / shape)
        )
    )
    expected = np.zeros(len(result.times))
    for step in range(1, len(expected)):
        expected[step] = initial * kernel[step - 1] + np.dot(
            expected[1:step], kernel[: step - 1][::-1]
        )
    assert result.incidence.sum() == pytest.approx(expected.sum(), rel=0.04)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"timestep": 0},
        {"duration": 0.015},
        {"initial_infections": 0},
        {"max_infections": 1},
        {"isolation": -1},
        {"tracing": 2},
    ],
)
def test_invalid_settings(kwargs):
    with pytest.raises(ValueError):
        simulate(1, 1, 2, **kwargs)
