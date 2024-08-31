import numpy as np
from scipy.integrate import quad as QUAD
import matplotlib.pyplot as plt
from scipy.stats import linregress
from timeit import timeit

from halton_sequences import generate_halton_vectors
from halton_sequences import halton_sequence
from halton_sequences import apply_box_muller_transform


def monte_carlo(num_paths, payoff_function, initial_stock_price, expiration_time, interest_rate, mean_function, deviation_function, method="normal", halton_primes=(None, None)):
    """
    Perform Monte Carlo simulation to estimate the option value.
    
    Parameters:
    - num_paths: Number of simulation paths
    - payoff_function: Function to calculate the option payoff
    - initial_stock_price: Initial stock price
    - expiration_time: Time to expiration (in years)
    - interest_rate: Risk-free interest rate
    - mean_function: Function to calculate the mean stock price
    - deviation_function: Function to calculate the standard deviation of the stock price
    - method: Method for generating random variables ('normal', 'antithetic', 'moment_matching', 'halton')
    - halton_primes: Tuple of Halton sequence primes (optional)
    
    Returns:
    - option_value: Estimated option value
    - option_error: Standard error of the option value
    """
    
    # Initialize parameters
    rng = np.random.default_rng(seed=0)
    
    if method == "normal":
        random_vars = rng.normal(0.0, 1.0, size=num_paths)
    
    elif method == "antithetic":
        random_vars = rng.normal(0.0, 1.0, size=num_paths // 2)
        random_vars = np.concatenate([random_vars, -random_vars])
    
    elif method == "moment_matching":
        random_vars = rng.normal(0.0, 1.0, size=num_paths // 2)
        variance = np.var(random_vars)
        normalized_vars = random_vars / np.sqrt(variance)
        random_vars = np.concatenate([normalized_vars, -normalized_vars])
    
    elif method == "halton":
        halton_prime_a, halton_prime_b = halton_primes
        halton_sequence = generate_halton_vectors(halton_prime_a, halton_prime_b, num_paths)
        random_vars = apply_box_muller_transform(halton_sequence)
    
    else:
        print("Error: Invalid method specified for random variable generation.")
        return None, None
    
    # Calculate terminal stock prices
    terminal_stock_prices = mean_function(initial_stock_price, expiration_time) + deviation_function(initial_stock_price, expiration_time) * np.sqrt(expiration_time) * random_vars
    
    # Calculate payoff values
    payoff_values = payoff_function(terminal_stock_prices)
    
    # Estimate option value and error
    option_value = np.mean(payoff_values) * np.exp(-interest_rate * expiration_time)
    option_error = np.sqrt(np.var(payoff_values)) / np.sqrt(num_paths)
    
    return option_value, option_error


def path_dependent_option(num_paths, num_steps, initial_stock_price, theta, alpha, beta, gamma, sigma, interest_rate, expiration_time):
    """
    Calculate the value of a path-dependent option using Monte Carlo simulation.
    
    Parameters:
    - num_paths: Number of simulation paths
    - num_steps: Number of time steps in the simulation
    - initial_stock_price: Initial stock price
    - theta: Drift term for the stock price
    - alpha: Parameter for drift term adjustment
    - beta: Parameter for drift term adjustment
    - gamma: Parameter for the volatility term
    - sigma: Base volatility of the stock
    - interest_rate: Risk-free interest rate
    - expiration_time: Time to expiration (in years)
    
    Returns:
    - option_value: Estimated option value
    - option_error: Standard error of the option value
    """
    # Initialize parameters
    rng = np.random.default_rng(seed=0)
    random_vars = rng.normal(0.0, 1.0, size=(num_paths, num_steps))
    
    # Adjust volatility if additional volatility is provided
    delta_time = expiration_time / num_steps
    
    # Initialize stock price paths
    stock_paths = np.zeros((num_paths, num_steps))
    stock_paths[:, 0] = initial_stock_price
    
    # Simulate stock price paths
    for step in range(num_steps - 1):
        stock_paths[:, step + 1] = (stock_paths[:, step] + 
                                    (alpha * theta - beta * stock_paths[:, step]) * delta_time + 
                                    sigma * np.abs(stock_paths[:, step]) ** gamma * np.sqrt(delta_time) * random_vars[:, step])
    
    # Calculate the option payoff
    min_stock_prices = np.min(stock_paths[:, 1:], axis=1)
    payoffs = np.maximum(0, stock_paths[:, -1] - min_stock_prices)
    
    # Estimate option value and error
    option_value = np.exp(-interest_rate * expiration_time) * np.mean(payoffs)
    option_error = np.sqrt(np.var(payoffs)) / np.sqrt(num_paths)
    
    return option_value, option_error

def path_dependent_derivative(num_paths, num_steps, initial_stock_price, theta, alpha, beta, gamma, sigma, interest_rate, expiration_time, delta_sigma):
    """
    Estimate the derivative of the option value with respect to the parameter sigma.
    
    Parameters:
    - num_paths: Number of simulation paths
    - num_steps: Number of time steps in the simulation
    - initial_stock_price: Initial stock price
    - theta: Drift term for the stock price
    - alpha: Parameter for drift term adjustment
    - beta: Parameter for drift term adjustment
    - gamma: Parameter for the volatility term
    - sigma: Base volatility of the stock
    - interest_rate: Risk-free interest rate
    - expiration_time: Time to expiration (in years)
    - delta_stock_price: Small change in initial stock price
    
    Returns:
    - option_derivative: Estimated derivative of the option value with respect to the initial stock price
    """
    # Calculate option value with initial stock price
    value_with_ds, _ = path_dependent_option(num_paths, num_steps, initial_stock_price, theta, alpha, beta, gamma, sigma, interest_rate, expiration_time)
    
    # Calculate option value without additional change
    value_without_ds, _ = path_dependent_option(num_paths, num_steps, initial_stock_price, theta, alpha, beta, gamma, sigma+delta_sigma, interest_rate, expiration_time)
    
    # Estimate derivative
    option_derivative = (value_with_ds - value_without_ds) / delta_sigma
    
    return option_derivative