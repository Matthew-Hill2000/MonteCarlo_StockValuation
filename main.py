import numpy as np
from scipy.integrate import quad as QUAD
import matplotlib.pyplot as plt
from scipy.stats import linregress
from timeit import timeit

from monte_carlo_methods import monte_carlo
from analysis_functions import *

plt.rcParams.update({'font.size': 15})

START_TIME = 0
INITIAL_STOCK_PRICE = 70154
INTEREST_RATE = 0.01
X_1 = 70000
X_2 = 80000
# MODEL PARAMETERS
THETA = 70100
ALPHA = 0.01
BETA = 0.01
GAMMA = 1.05
SIGMA = 0.26

def payoff_function(S):
    """
    Calculate the payoff of an option based on the underlying stock price.

    Parameters:
    - S: A numpy array or scalar representing the stock price(s).

    Returns:
    - A numpy array or scalar representing the option payoff for the given stock price(s).
    """
    return np.where(S < X_1, X_2 - S, 
                    np.where(S < X_2, S - X_2, 
                             X_1 - S))


def calculate_analytical_value(initial_stock_price, expiration_date, mean_function, variance_function, payoff_function, interest_rate):
    """
    Calculate the exact value of the option using analytical integration.

    Parameters:
    - initial_stock_price: The initial stock price (S0).
    - expiration_date: The expiration time of the option (T).
    - mean_function: A function that returns the mean of the stock price distribution.
    - variance_function: A function that returns the variance of the stock price distribution.
    - payoff_function: A function that computes the option payoff given the stock price.
    - interest_rate: The risk-free interest rate (r).

    Returns:
    - value: The calculated exact option value.
    - error: The estimated error of the integration.
    """
    S_0 = initial_stock_price
    T = expiration_date
    f = mean_function
    v = variance_function
    g = payoff_function
    r = interest_rate

    # Define the integrand function for numerical integration
    C_integrand = lambda z: g(z) * np.exp(-0.5 * (z - f(S_0, T))**2 / (v(S_0, T)**2 * T))
    
    # Perform the integration in three intervals
    I1, I2, I3 = QUAD(C_integrand, -9900000.0, X_1), QUAD(C_integrand, X_1, X_2), QUAD(C_integrand, X_2, 9900000.0)
    
    # Compute the exact option value and error
    value = np.exp(-r * T) * (I1[0] + I2[0] + I3[0]) / (v(S_0, T) * np.sqrt(2 * np.pi * T))
    error = I1[1] + I2[1] + I3[1]
    
    return value, error

def single_monte_carlo():
    """
    Main function to perform Monte Carlo simulation and compare it with the analytical value.
    
    This function sets the parameters for the simulation, runs the Monte Carlo method to estimate the option value,
    and calculates the exact value using the analytical approach. It then prints the results.
    """
    T = 2.0
    N = 1000000

    # Define the mean and variance functions
    f = lambda S_0, T: S_0 * (np.cosh(2 * BETA * T - ALPHA * T) - 1) + THETA * (3 - np.exp(ALPHA * T) - np.exp(BETA * T))
    v = lambda S_0, T: SIGMA * (1 + ALPHA * T) * 0.5 * (S_0 + THETA)**GAMMA

    # Perform Monte Carlo simulation
    value = monte_carlo(N, payoff_function, INITIAL_STOCK_PRICE, T, INTEREST_RATE, f, v, method="moment_matching", halton_primes=(2, 3))
    print("Monte Carlo Value:", value)

    # Calculate exact analytical value
    value_exact = calculate_analytical_value(INITIAL_STOCK_PRICE, T, f, v, payoff_function, INTEREST_RATE)
    print("Analytical Value:", value_exact)



def time_investigation():
    num_paths = 10000
    num_runs = 100
    T = 2.0

    # Define the mean and variance functions
    f = lambda S_0, T: S_0 * (np.cosh(2 * BETA * T - ALPHA * T) - 1) + THETA * (3 - np.exp(ALPHA * T) - np.exp(BETA * T))
    v = lambda S_0, T: SIGMA * (1 + ALPHA * T) * 0.5 * (S_0 + THETA)**GAMMA

    # Run the time analysis
    time_analysis(num_paths, payoff_function, INITIAL_STOCK_PRICE, T, INTEREST_RATE, f, v, halton_primes=(2,3), num_runs=num_runs)

def error_investigation():
    
    T = 2.0
    # Define the mean and variance functions
    f = lambda S_0, T: S_0 * (np.cosh(2 * BETA * T - ALPHA * T) - 1) + THETA * (3 - np.exp(ALPHA * T) - np.exp(BETA * T))
    v = lambda S_0, T: SIGMA * (1 + ALPHA * T) * 0.5 * (S_0 + THETA)**GAMMA

    analytical_value, _ = calculate_analytical_value(INITIAL_STOCK_PRICE, T, f, v, payoff_function, INTEREST_RATE)

    N_analysis(analytical_value, payoff_function, INITIAL_STOCK_PRICE, T, INTEREST_RATE, f, v, halton_primes=(2, 3))

    path_dep_N_analysis()

if __name__ == "__main__":
    single_monte_carlo()

    time_investigation()

    error_investigation()


