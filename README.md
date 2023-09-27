# StockValuation

This report concerns the calculation of various financial products using a variety of Monte Carlo methods. The first task was to calculate the value of a European put option using a standard Monte Carlo method, followed by the implementation of antithetic variables, moment matching and Halton sequences. The basis of the Monte Carlo method is to simulate many separate paths for the share price, from which we can calculate the option value in each case, taking the average to be our expected value. The simulations are run $N$ number of times, from which the average option value is calculated. This is then repeated $M$ number of times, forming a normally distributed set of Option values and averaged again to find the expected value. This is achieved by using techniques for drawing random variables $\phi_i$ from a normal distribution to simulate the statistical uncertainty in the evolution of the stock price. The initial Monte Carlo method simply uses the pseudo random number generator mt19937 to take draws from the standard normal distribution class of c++. \\

This is improved upon by utilizing the technique of antithetic variables. for every random draw $\phi_i$ from the normal distribution, we also take its negative. This ensures we have a distribution of random draws with mean of zero and hence the mean of the sample paths is also correct. We therefore only have to take $N/2$ random draws to result in $N$ Values. \\

Further to this we can then use moment matching to ensure that the the variance of the sample paths match that of the required distribution by forcing the variance of our random draws from the normal distribution to equal 1. We start by taking our random draws $\phi_i$ using the antithetic technique, and calculate the unbiased estimator for the variance $\nu^2$


$$\frac{2 \sum_i \phi_i^2}{2n-1} = \nu^2$$
