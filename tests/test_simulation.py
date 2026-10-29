import numpy as np
import pytest
from scipy.stats import gamma

import fraser.simulation as simulation_module
from fraser.simulation import (
    _People,
    _process_removals,
    binned_incidence,
    estimate_growth,
    find_critical_r0,
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


def test_exponential_growth_fit():
    times = np.arange(0.5, 10)
    assert estimate_growth(times, 200 * np.exp(0.3 * times)) == pytest.approx(0.3)


def test_growth_fit_includes_first_zero_as_one_and_ignores_later_bins():
    times = np.arange(0.5, 4)
    incidence = np.array([100, 10, 0, 1000])
    assert estimate_growth(times, incidence) == pytest.approx(-np.log(10))
    assert incidence[2] == 0


def test_growth_fit_with_only_one_bin_is_unavailable():
    assert np.isnan(estimate_growth(np.array([0.5]), np.array([0])))


def test_incidence_binning_excludes_initial_seeds():
    result = simulate(0, 1, 2, duration=2)
    times, incidence = binned_incidence(result)
    np.testing.assert_array_equal(times, [0.5, 1.5])
    np.testing.assert_array_equal(incidence, [0, 0])


def test_partial_incidence_bin_preserves_counts():
    result = simulate(1.5, 1, 2, initial_infections=100, duration=2.5)
    times, incidence = binned_incidence(result)
    np.testing.assert_array_equal(times, [0.5, 1.5, 2.25])
    assert np.dot(incidence, [1, 1, 0.5]) == result.incidence.sum()


def _mock_growth(monkeypatch, root, noise=(0, 0, 0)):
    calls = []

    def fake_simulate(r0, *args, **kwargs):
        calls.append(kwargs["seed"])
        return r0 - root + noise[(len(calls) - 1) % len(noise)]

    monkeypatch.setattr(simulation_module, "simulate", fake_simulate)
    monkeypatch.setattr(
        simulation_module,
        "binned_incidence",
        lambda rate, width: (np.array([0.0, 1.0]), np.exp([0.0, rate])),
    )
    return calls


def test_critical_bisection_brackets_zero(monkeypatch):
    calls = _mock_growth(monkeypatch, 3.7)
    result = find_critical_r0(1, 2)
    assert result.status == "bracketed"
    assert result.lower <= 3.7 <= result.upper
    assert result.upper - result.lower <= 0.25
    assert len(set(calls)) == len(calls)
    assert all(len(rates) == 3 for _, rates in result.evaluations)


def test_critical_median_tolerates_one_outlier(monkeypatch):
    _mock_growth(monkeypatch, 3.7, noise=(0, 0, 10))
    result = find_critical_r0(1, 2)
    assert result.lower <= 3.7 <= result.upper


@pytest.mark.parametrize("root,status", [(0, "below_range"), (20, "above_range")])
def test_critical_outside_bracket(monkeypatch, root, status):
    _mock_growth(monkeypatch, root)
    result = find_critical_r0(1, 2)
    assert result.status == status
    assert np.isnan(result.estimate)


def test_critical_unavailable_fit_is_unresolved(monkeypatch):
    _mock_growth(monkeypatch, 3.7)
    monkeypatch.setattr(simulation_module, "estimate_growth", lambda *args: np.nan)
    assert find_critical_r0(1, 2).status == "unresolved"
