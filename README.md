# Theory
## Stock Options

When the price of a stock is said to follow a risk-neutral distribution, we have that its price at a time $t$ is given by

$$ S_t \approx N(f(S_0,t), v^2(S_0, t)t) $$

where $S_0$ is the current stock price and $f(S_0,t)$ and $v(S_0,t)$ are functions calibrated for the specific context. For the financial contract to be valued in this report, the price at time $t$ is given by

$$S_t = S_0 (\cosh(2 \beta T - \alpha T)-1) + \theta (3-e^{\alpha T} - e^{\beta T}) + \sigma (1 + \alpha T) \frac{1}{2} (S_0 + \theta)^{\gamma} \sqrt{T} \phi$$

where $\theta=70100$, $\alpha=0.01$, $\beta=0.01$, $\gamma=1.05$, $\sigma=0.26$, $T=2$ and  $\phi = N(0,1)$.

If we consider a financial contract $C(S,T)$ written on the underlying stock $S$, which has a payoff given by

$$
    C(S,T) = g(S) = 
    \begin{cases} 
X_2 - S_T & \text{if } S_T < X_1 \\
S_T - X_2 & \text{if } X_1 \leq S_T < X_2 \\
X_1 - S_T & \text{if } S_T \geq X_2 
\end{cases}
$$

Then the analytical solution for any given payoff is given by the numerical integration

$$C(S_0, t=0) = \frac{e^{-rT}}{v \sqrt{2 \pi T}} \int_{-\infty}^{\infty} g(z) \exp  \left[ -\frac{(z-f)^2}{2v^2T} \right]  dz$$

We can carry out a Monte-Carlo valuation for the option by sampling from the normal distribution $\phi$ and calculating the stock price several times, with the $i^{th}$ stock price given by

$$S_T^i = f(S_0, T) + v(S_0, T) \sqrt{T} \phi_i$$

We then average over the $n$ calculated stock prices to approximate the value of the financial contract as

$$C(S_0, t=0) \approx e^{-rT} \frac{1}{n} \sum_{i=1}^{n} g(S_T^i)$$

## Path Dependent Options

Assuming that the risk-neutral stochastic process follows the SDE

$$ds = f(S,t)dt + v(S,t)dW.$$

A path-dependent option depends on all of the share prices, $S(t_k)$, at $K+1$ equally spaced sampling times $t_0, t_1,...,t_k=k \Delta t,..., t_K$ with $t_0=0$, $t_K=T$ and

$$\Delta t = \frac{T}{K}.$$

In this report, the SDE in equation (8) takes the specific form

$$ds = (\alpha \theta - \beta S)dt + \sigma (|S|)^{\gamma}dW,$$

where $\theta=70100$, $\alpha=0.01$, $\beta=0.01$, $\gamma=1.05$, $\sigma=0.26$, $T=1$.

Since the options are path-dependent, We must approximate the price at each of K time steps between the initial time and the time of maturity $T$. If the time step is given by $\Delta T$, then the share price at each discreet point in time is given by

$$S^i(t_k) = S^i(t_{k-1}) + f(S^i(t_{k-1}), t_{k-1}) \Delta t + v(S^i(t_{k-1}), t_{k-1})\sqrt{\Delta t} \phi_{i,k-1}$$

In this report, we value a minimum floating strike lookback call option, the payoff, $G$, is found by first calculating the minimum share price

$$A = \min_k S(t_k)$$

and then 

$$G(S,A) = \max(S-A, 0)$$

To calculate the value of the contract at $t=0$ using Monte Carlo simulation, we average over $n$ approximations and apply a discounting factor to get

$$C(S_0, t=0) \approx e^{-rT} \frac{1}{n} \sum_{i=1}^{n} G(S^i, A)$$







