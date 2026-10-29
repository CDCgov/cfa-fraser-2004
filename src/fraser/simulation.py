"""One-population timestep branching process with symptom-triggered tracing."""

import math
from collections import deque
from dataclasses import dataclass

import numpy as np
from scipy.special import betainc
from scipy.stats import gamma


def _positive(value: float, name: str) -> None:
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")


def presymptomatic_fraction(symptom_mean: float, shape: float) -> float:
    """Derive theta for common-shape gammas with mean infectiousness time 1."""
    _positive(symptom_mean, "symptom_mean")
    _positive(shape, "shape")
    return float(betainc(shape, shape, symptom_mean / (1 + symptom_mean)))


@dataclass
class SimulationResult:
    """Per-step new events; incidence excludes the initial infections.

    A population-limit result contains only fully completed timesteps.
    Isolation and quarantine count events, so one person can appear in both.
    """

    times: np.ndarray
    incidence: np.ndarray
    isolations: np.ndarray
    quarantines: np.ndarray
    r0: float
    symptom_mean: float
    shape: float
    isolation: float
    tracing: float
    theta: float
    initial_infections: int
    total_infections: int
    timestep: float
    requested_duration: float
    seed: int | None
    completed: bool


class _People:
    def __init__(self, capacity: int):
        self.size = 0
        self.infection_step = np.empty(capacity, dtype=np.int64)
        self.symptom_step = np.empty(capacity, dtype=np.int64)
        self.parent = np.empty(capacity, dtype=np.int64)
        self.first_child = np.empty(capacity, dtype=np.int64)
        self.next_sibling = np.empty(capacity, dtype=np.int64)
        self.isolated = np.empty(capacity, dtype=bool)
        self.quarantined = np.empty(capacity, dtype=bool)
        self.symptoms_processed = np.empty(capacity, dtype=bool)
        self.contacts_traced = np.empty(capacity, dtype=bool)

    def add(self, parents: np.ndarray, step: int, symptoms: np.ndarray) -> None:
        count = len(parents)
        if not count:
            return
        end = self.size + count
        capacity = len(self.parent)
        if end > capacity:
            for name, old in list(vars(self).items()):
                if isinstance(old, np.ndarray):
                    enlarged = np.empty(max(end, 2 * capacity), dtype=old.dtype)
                    enlarged[: self.size] = old[: self.size]
                    setattr(self, name, enlarged)
        newborns = slice(self.size, end)
        self.infection_step[newborns] = step
        self.symptom_step[newborns] = symptoms
        self.parent[newborns] = parents
        self.first_child[newborns] = -1
        self.next_sibling[newborns] = -1
        self.isolated[newborns] = False
        self.quarantined[newborns] = False
        self.symptoms_processed[newborns] = False
        self.contacts_traced[newborns] = False
        if parents[0] >= 0:
            starts = np.r_[0, np.flatnonzero(np.diff(parents)) + 1]
            ends = np.r_[starts[1:] - 1, count - 1]
            parent_ids = parents[starts]
            self.next_sibling[newborns] = np.arange(self.size + 1, end + 1)
            self.next_sibling[self.size + ends] = self.first_child[parent_ids]
            self.first_child[parent_ids] = self.size + starts
        self.size = end


def _process_removals(
    people: _People,
    step: int,
    isolation: float,
    tracing: float,
    rng: np.random.Generator,
) -> tuple[int, int]:
    current = slice(0, people.size)
    due = np.flatnonzero(
        ~people.symptoms_processed[current] & (people.symptom_step[current] <= step)
    )
    people.symptoms_processed[due] = True
    selected = due[people.quarantined[due] | (rng.random(len(due)) < isolation)]
    people.isolated[selected] = True
    new_isolations = len(selected)
    new_quarantines = 0
    if tracing == 0:
        return new_isolations, new_quarantines
    queue = deque(selected.tolist())
    while queue:
        person = queue.popleft()
        if people.contacts_traced[person]:
            continue
        people.contacts_traced[person] = True
        children = []
        child = people.first_child[person]
        while child >= 0:
            if not (people.isolated[child] or people.quarantined[child]):
                children.append(child)
            child = people.next_sibling[child]
        contacts = np.asarray(children, dtype=np.int64)
        contacts = contacts[rng.random(len(contacts)) < tracing]
        people.quarantined[contacts] = True
        new_quarantines += len(contacts)
        symptomatic = contacts[people.symptoms_processed[contacts]]
        people.isolated[symptomatic] = True
        new_isolations += len(symptomatic)
        queue.extend(symptomatic.tolist())
    return new_isolations, new_quarantines


def simulate(
    r0: float,
    symptom_mean: float,
    shape: float,
    *,
    isolation: float = 1.0,
    tracing: float = 1.0,
    initial_infections: int = 1_000,
    timestep: float = 0.01,
    duration: float = 10.0,
    max_infections: int = 1_000_000,
    seed: int | None = 2004,
) -> SimulationResult:
    """Simulate many independent seed trees together, with mean generation time 1.

    Transmissions precede symptom processing in every step. Quarantine blocks
    transmission but retains symptoms and subsequent tracing. The horizon must
    be an integer number of timesteps. No growth rate or threshold is estimated.
    """
    theta = presymptomatic_fraction(symptom_mean, shape)
    if not math.isfinite(r0) or r0 < 0:
        raise ValueError("r0 must be nonnegative and finite")
    for value, name in [(isolation, "isolation"), (tracing, "tracing")]:
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must be between 0 and 1")
    _positive(timestep, "timestep")
    _positive(duration, "duration")
    steps = round(duration / timestep)
    if steps < 1 or not math.isclose(steps * timestep, duration):
        raise ValueError("duration must be an integer number of timesteps")
    if not isinstance(initial_infections, int) or initial_infections < 1:
        raise ValueError("initial_infections must be a positive integer")
    if not isinstance(max_infections, int) or max_infections < initial_infections:
        raise ValueError("max_infections must be an integer >= initial_infections")
    rng = np.random.default_rng(seed)
    people = _People(max(initial_infections, 16))

    def add_people(parents: np.ndarray, step: int) -> None:
        delays = rng.gamma(shape, scale=symptom_mean / shape, size=len(parents))
        symptom_steps = step + np.maximum(1, np.ceil(delays / timestep)).astype(
            np.int64
        )
        people.add(parents, step, symptom_steps)

    add_people(np.full(initial_infections, -1, dtype=np.int64), 0)
    means = (
        r0
        * gamma.pdf(np.arange(steps + 1) * timestep, shape, scale=1 / shape)
        * timestep
    )
    incidence = np.zeros(steps + 1, dtype=np.int64)
    isolations = np.zeros(steps + 1, dtype=np.int64)
    quarantines = np.zeros(steps + 1, dtype=np.int64)
    completed_steps = 0
    for step in range(1, steps + 1):
        current = slice(0, people.size)
        eligible = np.flatnonzero(
            ~people.isolated[current] & ~people.quarantined[current]
        )
        counts = rng.poisson(means[step - people.infection_step[eligible]])
        births = int(counts.sum())
        if people.size + births > max_infections:
            break
        parents = np.repeat(eligible, counts)
        add_people(parents, step)
        incidence[step] = births
        isolations[step], quarantines[step] = _process_removals(
            people, step, isolation, tracing, rng
        )
        completed_steps = step
    recorded = slice(0, completed_steps + 1)
    return SimulationResult(
        times=np.arange(completed_steps + 1) * timestep,
        incidence=incidence[recorded],
        isolations=isolations[recorded],
        quarantines=quarantines[recorded],
        r0=r0,
        symptom_mean=symptom_mean,
        shape=shape,
        isolation=isolation,
        tracing=tracing,
        theta=theta,
        initial_infections=initial_infections,
        total_infections=people.size,
        timestep=timestep,
        requested_duration=duration,
        seed=seed,
        completed=completed_steps == steps,
    )
