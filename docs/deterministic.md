# Deterministic mathematics

The paper begins with definitions that culminate in Equation 6, which describes the [isolation-only](#isolation-only) case.
A general derivation of contact tracing culminates in Equation 16.
Strong simplifying assumptions produce the [tracing and quarantine](#tracing-and-quarantining) result in Equation 18.

## Isolation only

Let $\tau$ denote time since infection.
The infectiousness density is $\beta(\tau) = R_0 g(\tau)$, with normalized density $g(\tau)$ such that $\int_0^\infty g(\tau)\,d\tau = 1$.
The mean (no-intervention counterfactual) generation interval is $T_g = \int_0^\infty\tau g(\tau)\,d\tau$.

Let $S(\tau)$ be the incubation survival function (i.e., the probability symptoms have not appeared by $\tau$).
(The asymptomatic fraction is $S(\infty)$; this is zero for the exponential example below.)
Then the fraction $\theta$ of transmission before symptom onset (including asymptomatic infections) is:

$$
  \theta = \int_0^\infty g(\tau)S(\tau)\,d\tau
$$

Isolation leaves presymptomatic transmission intact and reduces symptomatic transmission by $\varepsilon_I$.

The key relation is Equation 6, which gives the effective reproduction number with isolation:

$$
  R_I = R_0(1 - \varepsilon_I + \varepsilon_I\theta)
$$

Control requires $R_I < 1$.
With perfect isolation $\varepsilon_I = 1$, this is $\theta < 1/R_0$.

## Tracing and quarantining

Equation 16 gives a general result, when assuming that isolation and contact tracing are independent events.
Further assuming that infectiousness and symptom onset are exponential-distributed gives the highly tractable Equation 18 for the control threshold:

$$
  1
  = R_0 \left\{ \left[ (1 - \varepsilon_I)(1 - \varepsilon_T) + \varepsilon_I(1 - \varepsilon_T)\theta + (1 - \varepsilon_I)\varepsilon_T \theta + \varepsilon_I\varepsilon_T \frac{\theta}{2 - \theta} \right] \right\}
$$

where $\varepsilon_T$ is the efficacy of contact tracing.

### Equation 18 typo

Equation 18 prints $2 + \theta$; it should be $2 - \theta$.
At $\theta = 1$ the corrected formula gives threshold $R_0 = 1$ regardless of efficacy; the printed formula does not.

Equation 17 gives the threshold condition:

$$
  \begin{align*}
    1 & = R_0\int_0^\infty g(t) [1 - \varepsilon_I + \varepsilon_I S(t)] [1 - \varepsilon_T + \varepsilon_T S(t)]\,dt                                                                                                            \\
      & = R_0\left[ (1 - \varepsilon_I)(1 - \varepsilon_T) + \{\varepsilon_I(1 - \varepsilon_T) + (1 - \varepsilon_I)\varepsilon_T\} \int_0^\infty g(t)S(t)\,dt + \varepsilon_I\varepsilon_T\int_0^\infty g(t)S(t)^2\,dt \right]
  \end{align*}
$$

Assuming exponential incubation times $S(t) = e^{-\nu t}$ and exponential infectiousness $g(t) = e^{-t}$ (choosing time units so that $T_g = 1$), then the first integral is simply

$$
  \theta = \int_0^\infty g(t)S(t)\,dt
         = \int_0^\infty e^{-(1+\nu)t}\,dt = \frac{1}{1 + \nu}.
$$

The second is

$$
  \int_0^\infty g(t)S(t)^2\,dt = \int_0^\infty e^{-(1+2\nu)t}\,dt
                               = \frac{1}{1 + 2\nu} = \frac{\theta}{2 - \theta}.
$$
