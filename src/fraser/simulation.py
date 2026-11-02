"""One-population timestep branching process with symptom-triggered tracing."""

import math
from collections import deque
from dataclasses import dataclass, field

import numpy as np
import numpy.random
from scipy.special import betainc
from scipy.stats import gamma


def presymptomatic_fraction(symptom_mean: float, shape: float) -> float:
    """Derive theta for common-shape gammas with mean infectiousness time 1."""
    return float(betainc(shape, shape, symptom_mean / (1 + symptom_mean)))


@dataclass(kw_only=True)
class CommonParameters:
    symptom_mean: float
    shape: float
    isolation: float = 0.9
    tracing: float = 1
    initial_infections: int = 1_000
    timestep: float = 0.01
    duration: float = 10
    max_infections: int = 1_000_000
    bin_width: float = 1
    seed: int = 2004
    # derived parameters
    theta: float = field(init=False)
    steps: int = field(init=False)

    def __post_init__(self):
        """Validate and derive parameters"""
        if not self.symptom_mean > 0:
            raise ValueError("symptom_mean must be positive and finite")
        if not self.shape > 0:
            raise ValueError("shape must be positive and finite")
        if not self.timestep > 0:
            raise ValueError("timestep must be positive and finite")
        if not self.duration > 0:
            raise ValueError("duration must be positive and finite")
        if not 0 <= self.isolation <= 1:
            raise ValueError("isolation must be between 0 and 1")
        if not 0 <= self.tracing <= 1:
            raise ValueError("tracing must be between 0 and 1")
        if self.initial_infections <= 0:
            raise ValueError("initial_infections must be positive")
        if self.timestep <= 0:
            raise ValueError("timestep must be positive")
        if self.duration <= 0:
            raise ValueError("duration must be positive")
        if self.max_infections <= 0:
            raise ValueError("max_infections must be positive")
        if not isinstance(self.initial_infections, int) or self.initial_infections < 1:
            raise ValueError("initial_infections must be a positive integer")
        if (
            not isinstance(self.max_infections, int)
            or self.max_infections < self.initial_infections
        ):
            raise ValueError("max_infections must be an integer >= initial_infections")
        if self.bin_width <= 0:
            raise ValueError("bin_width must be positive")

        self.steps = round(self.duration / self.timestep)
        self.theta = presymptomatic_fraction(self.symptom_mean, self.shape)

        if self.steps < 1 or not math.isclose(
            self.steps * self.timestep, self.duration
        ):
            raise ValueError("duration must be an integer number of timesteps")


@dataclass(kw_only=True)
class SimulationParameters(CommonParameters):
    r0: float

    def __post_init__(self):
        """Validate R0 and the shared simulation parameters."""
        if not math.isfinite(self.r0) or self.r0 < 0:
            raise ValueError("r0 must be nonnegative and finite")
        super().__post_init__()


@dataclass(kw_only=True)
class SearchParameters(CommonParameters):
    lower: float = 1
    upper: float = 15
    tolerance: float = 0.25
    repeats: int = 3

    def __post_init__(self):
        """Validate the search parameters and shared configuration."""
        super().__post_init__()
        if self.lower >= self.upper:
            raise ValueError("require lower < upper")
        if self.tolerance <= 0 or self.repeats < 1:
            raise ValueError("require positive tolerance and repeats >= 1")

    def simulation_parameters(self, r0: float, seed: int) -> SimulationParameters:
        """Build one simulation candidate for this search."""
        return SimulationParameters(
            r0=r0,
            symptom_mean=self.symptom_mean,
            shape=self.shape,
            isolation=self.isolation,
            tracing=self.tracing,
            initial_infections=self.initial_infections,
            timestep=self.timestep,
            duration=self.duration,
            max_infections=self.max_infections,
            bin_width=self.bin_width,
            seed=seed,
        )


@dataclass
class SimulationResult:
    """Per-step new events; incidence excludes the initial infections.

    A population-limit result contains only fully completed timesteps.
    """

    parameters: SimulationParameters
    times: np.ndarray
    incidence: np.ndarray
    completed: bool

    def binned_incidence(
        self,
        bin_width: float = 1.0,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return bin midpoints and new-infection rates, excluding seed infections."""
        # number of timesteps per bin
        per_bin = round(bin_width / self.parameters.timestep)
        if per_bin < 1 or not math.isclose(
            per_bin * self.parameters.timestep, bin_width
        ):
            raise ValueError("bin_width must be an integer number of timesteps")

        # start and end times for each bin
        starts = np.arange(0, len(self.incidence) - 1, per_bin)
        ends = np.minimum(starts + per_bin, len(self.incidence) - 1)
        midpoints = (self.times[starts] + self.times[ends]) / 2
        if not len(starts):
            return np.array([]), np.array([])

        # Return rates so a final partial bin contributes proportionally when
        # callers integrate the binned series.
        counts = np.add.reduceat(self.incidence[1:], starts)
        widths = self.times[ends] - self.times[starts]
        return midpoints, counts / widths


def estimate_growth(times: np.ndarray, incidence: np.ndarray) -> float:
    """Fit log incidence from the start through the first zero, replacing it by 1.

    Returns NaN when fewer than two bins are available. This is an exploratory
    slope, not an uncertainty-adjusted or asymptotic growth estimate.
    """
    values = np.asarray(incidence, dtype=float).copy()
    zeros = np.flatnonzero(values == 0)
    if len(zeros):
        values = values[: zeros[0] + 1]
        values[-1] = 1
    if len(values) < 2:
        return math.nan
    return float(np.polyfit(times[: len(values)], np.log(values), 1)[0])


@dataclass
class CriticalResult:
    """A rough critical-R0 bracket, not a statistical confidence interval."""

    parameters: SearchParameters
    lower: float
    upper: float
    status: str

    @property
    def estimate(self) -> float:
        return (self.lower + self.upper) / 2 if self.status == "bracketed" else math.nan


def find_critical_r0(
    parameters: SearchParameters,
) -> CriticalResult:
    """Bisect median rough growth across independent seeds at each candidate R0.

    Results outside the initial bracket or with unavailable fits are explicit.
    Finite fits from population-limited histories are used as in the rough grid.
    """
    if not isinstance(parameters, SearchParameters):
        raise TypeError("parameters must be a SearchParameters instance")

    seed_sequence = np.random.SeedSequence(parameters.seed)

    def growth(r0: float) -> float:
        rates = []
        for stream in seed_sequence.spawn(parameters.repeats):
            seed = stream.generate_state(1)[0]
            result = simulate(parameters.simulation_parameters(r0, int(seed)))
            times, incidence = result.binned_incidence()
            rates.append(estimate_growth(times, incidence))
        return float(np.median(rates))

    def finish(status: str) -> CriticalResult:
        return CriticalResult(
            parameters=parameters,
            lower=lower,
            upper=upper,
            status=status,
        )

    lower, upper = parameters.lower, parameters.upper
    low_rate, high_rate = growth(lower), growth(upper)
    if not (math.isfinite(low_rate) and math.isfinite(high_rate)):
        return finish("unresolved")
    if low_rate > 0 and high_rate < 0:
        return finish("unresolved")
    if low_rate > 0:
        return finish("below_range")
    if high_rate < 0:
        return finish("above_range")
    while upper - lower > parameters.tolerance:
        midpoint = (lower + upper) / 2
        rate = growth(midpoint)
        if not math.isfinite(rate):
            return finish("unresolved")
        if rate > 0:
            upper = midpoint
        else:
            lower = midpoint
    return finish("bracketed")


class People:
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
        """
        Args:
            parents: indices of parents for each new infection
            step: current simulation step
            symptoms: symptom onset steps for each new infection
        """
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
    people: People,
    step: int,
    isolation: float,
    tracing: float,
    rng: np.random.Generator,
) -> None:
    current = slice(0, people.size)
    due = np.flatnonzero(
        ~people.symptoms_processed[current] & (people.symptom_step[current] <= step)
    )
    people.symptoms_processed[due] = True
    selected = due[people.quarantined[due] | (rng.random(len(due)) < isolation)]
    people.isolated[selected] = True

    if tracing == 0:
        return
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
        symptomatic = contacts[people.symptoms_processed[contacts]]
        people.isolated[symptomatic] = True
        queue.extend(symptomatic.tolist())


def simulate(p: SimulationParameters) -> SimulationResult:
    """Simulate many independent seed trees together, with mean generation time 1.

    Transmissions precede symptom processing in every step. Quarantine blocks
    transmission but retains symptoms and subsequent tracing. The horizon must
    be an integer number of timesteps. No growth rate or threshold is estimated.
    """

    rng = np.random.default_rng(p.seed)
    people = People(p.initial_infections)

    def add_people(people: People, parents: np.ndarray, step: int) -> None:
        """Add new people, drawing symptom onset steps"""
        delays = rng.gamma(p.shape, scale=p.symptom_mean / p.shape, size=len(parents))
        symptom_steps = step + np.maximum(1, np.ceil(delays / p.timestep)).astype(
            np.int64
        )
        people.add(parents, step, symptom_steps)

    # initialize infections
    add_people(people, np.full(p.initial_infections, -1, dtype=np.int64), 0)

    # initialize the readout
    incidence = np.zeros(p.steps + 1, dtype=np.int64)

    # pre-compute mean number of secondary infections, per timestep since infector's infection
    means = (
        p.r0
        * gamma.pdf(np.arange(p.steps + 1) * p.timestep, p.shape, scale=1 / p.shape)
        * p.timestep
    )

    # begin the loop
    completed_steps = 0
    for step in range(1, p.steps + 1):
        # current size of the simulation
        current = slice(0, people.size)

        # effectively infective people
        eligible = np.flatnonzero(
            ~people.isolated[current] & ~people.quarantined[current]
        )
        # number of secondary infections for each eligible infector
        counts = rng.poisson(means[step - people.infection_step[eligible]])
        # total number of new infections
        births = int(counts.sum())
        if people.size + births > p.max_infections:
            break
        parents = np.repeat(eligible, counts)
        add_people(people, parents, step)
        incidence[step] = births
        _process_removals(people, step, p.isolation, p.tracing, rng)
        completed_steps = step
    recorded = slice(0, completed_steps + 1)

    return SimulationResult(
        parameters=p,
        times=np.arange(completed_steps + 1) * p.timestep,
        incidence=incidence[recorded],
        completed=completed_steps == p.steps,
    )
