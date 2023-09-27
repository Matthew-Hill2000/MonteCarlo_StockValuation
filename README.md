# StockValuation

This report concerns the calculation of various financial products using a variety of Monte Carlo methods. The first task was to calculate the value of a European put option using a standard Monte Carlo method, followed by the implementation of antithetic variables, moment matching and Halton sequences. The basis of the Monte Carlo method is to simulate many separate paths for the share price, from which we can calculate the option value in each case, taking the average to be our expected value. The simulations are run $N$ number of times, from which the average option value is calculated. This is then repeated $M$ number of times, forming a normally distributed set of Option values and averaged again to find the expected value. This is achieved by using techniques for drawing random variables $\phi_i$ from a normal distribution to simulate the statistical uncertainty in the evolution of the stock price. The initial Monte Carlo method simply uses the pseudo random number generator mt19937 to take draws from the standard normal distribution class of c++. 

This is improved upon by utilizing the technique of antithetic variables. for every random draw $\phi_i$ from the normal distribution, we also take its negative. This ensures we have a distribution of random draws with mean of zero and hence the mean of the sample paths is also correct. We therefore only have to take $N/2$ random draws to result in $N$ Values. 

Further to this we can then use moment matching to ensure that the the variance of the sample paths match that of the required distribution by forcing the variance of our random draws from the normal distribution to equal 1. We start by taking our random draws $\phi_i$ using the antithetic technique, and calculate the unbiased estimator for the variance $\nu^2$


$$\frac{2 \sum_i \phi_i^2}{2n-1} = \nu^2.$$


We then divide each $\phi_i$ by $\nu$ and the distribution of random draws is then equal to 1. Finally, a version of the Monte Carlo method was implemented where the Halton sequence was used to generate random numbers. This works by first taking two prime numbers $a$ and $b$ before representing a series of integers in the base of $a$ and $b$, which can be expressed as


$$i = \sum_{j=1}^m x_jy^{j-1}$$

where i is the integer, y is the base and $0 \le x_j < y$. The Halton numbers are then given by

$$h(i;y) = \sum_{j=1}^m x_jy^{-j}.$$

For Halton numbers produced from two different primes, we can pair up the Halton numbers from the same integer to produce uniformly distributed coordinates in the unit square. if we label these coordinates $(x_1, x_2)$, redefining the $x$ symbol, then we can us the Box-Muller method produce two normally distributed variables $y_1$ and $y_2$ from 


$$y_1 = \cos(2 \pi x_2) \sqrt{-2 \log(x_1)}, \hspace{0.5cm} y_2 = \sin(2 \pi x_1) \sqrt{-2 \log(x_2)}.$$

For the methods using anithetic variables, moment matching and Halton sequences we take only N/2 random draws, such that the final result is N sample paths, allowing us to compare to the standard Monte Carlo method. For the Halton sequence, i included a function to choose two random primes for every sample path such that the paths differed. The algorithm for the Halton sequence is outlined in psuedo code in \cite{11}.


#Stock Options
 
The value of the stock price at time t is claimed to be given by the risk-neutral distribution \cite{4}


$$S_t = N(f(S_0,t), v^2(S_0, t)t)$$

for some functions $f(S_0,t)$ and $v(S_0,t)$, which in this report is given by


$$S_t = S_0 (\alpha T + \tan(\beta T)) + \theta \cosh(2 \beta T - \alpha T) + \sigma (1 + \beta T)^2 S_0 e^{2 \gamma - 2 + \alpha T} \sqrt{T} \phi$$



where 

$$\phi = N(0,1).$$

The initial stock price is $S_0=500.72$, and the risk-free interest rate is $r=0.03$. The option matures at $T=2$ with a strike price of $X=500$. The parameters are $\theta=501$, $\alpha=0.01$, $\beta=0.01$, $\gamma=1.01$, and the volatility of the option is $\sigma=0.23$. The analytical solution to the value of the put option is given by


$$P(S_0, t=0) = [XN(z) + v(S_0,T) \sqrt{T} \frac{1}{\sqrt{2\pi}}e^{-z^2/2}-f(S_0,T)N(z)]e^{-rT}$$

where N(z) is the cumulative normal distribution and z is given by


$$z = \frac{X-f(S_0,T)}{v(S_0,T)\sqrt{T}}.$$

To calculate the valuation using the Monte Carlo method we take samples from the standard normal distribution to evaluate


$$S_T^i = f(S_0, T) + v(S_0, T) \sqrt{T} \phi_i$$


which then allows us to calculate the option value from

$$P(S_0, t=0) = e^{-rT} \frac{1}{n} \sum_{i=1}^{n} max(X-S_T^i, 0)$$

#Path Dependent Options

Assuming that the risk neutral stochastic process follows the SDE

$$ds = f(S,t)dt + v(S,t)dW.$$

For this report we take


$$ds = (\alpha \theta - \beta S)dt + \sigma (|S|)^{\gamma}dW,$$

where the parameters are the same as those from section 1.2, however now with T=1. With path dependent options the payout is dependent on $S(t_k)$, the share prices at K+1 equally spaced times from $t_0$ to $t_K$. The time step is given by

$$\Delta t = \frac{T}{K}.$$

Since the final stock price is not computed from a single jump, its distributions will not be normal. To calculate S for each subsequent time we use the equation


$$S^i(t_k) = S^i(t_{k-1}) + f(S^i(t_{k-1}), t_{k-1}) \Delta t + v(S^i(t_{k-1}), t_{k-1})\sqrt{\Delta t} \phi_{i,k-1}$$

In this report we value an Asian option, the payoff of which is given by first taking the average

$$A^i = \frac{1}{K} \sum_{k=1}^{K} S^i(t_k)$$

and then specifically for a fixed strike put option calculating the value from

$$P(S_0, t=0) = e^{-rT} \frac{1}{n} \sum_{i=1}^{n} max(X-A^i, 0)$$

where X is the strike price.
