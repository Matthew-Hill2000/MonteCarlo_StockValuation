import numpy as np
from scipy.integrate import quad as QUAD
import matplotlib.pyplot as plt
from scipy.stats import linregress
from timeit import timeit
from monte_carlo_methods import monte_carlo
from monte_carlo_methods import path_dependent_option
from monte_carlo_methods import path_dependent_derivative
from halton_sequences import halton_sequence
from halton_sequences import generate_halton_vectors
from halton_sequences import apply_box_muller_transform

 # Error Analysis
def error_analysis(N_values, errors):
    """
    Performs error analysis by plotting the Monte Carlo errors against N^{-1/2} 
    and fitting a linear regression line to the data.

    Parameters:
    - N_values: Array of sample sizes (N).
    - errors: Array of errors corresponding to the sample sizes.

    Returns:
    - A plot showing the Monte Carlo errors and the fitted linear regression line.
    """
    N_inv_sqrt = N_values**-0.5
    
    slope, intercept, _, _, _ = linregress(N_inv_sqrt, errors)
    fitted_errors = slope * N_inv_sqrt + intercept
    
    plt.figure(figsize=(10, 6))
    plt.plot(N_inv_sqrt, errors, 'o', label="Monte Carlo Errors")
    plt.plot(N_inv_sqrt, fitted_errors, 'r', label=f"Fitted Line (y = {slope:.2f}x + {intercept:.2f})")
    
    plt.xlabel("$N^{-1/2}$")
    plt.ylabel("Error")
    plt.title("Monte Carlo Error vs. $N^{-1/2}$")
    plt.legend()
    plt.grid(True)
    plt.show()

# Multi Error Analysis
def multi_error_analysis(N_values, monte_errors, anti_errors, moment_errors, halton_errors):
    """
    Performs error analysis for multiple Monte Carlo methods and plots the errors 
    against N^{-1/2} for comparison.

    Parameters:
    - N_values: Array of sample sizes (N).
    - monte_errors: Array of errors from the standard Monte Carlo method.
    - anti_errors: Array of errors from the antithetic variates method.
    - moment_errors: Array of errors from the moment matching method.
    - halton_errors: Array of errors from the Halton sequence method.

    Returns:
    - A plot comparing the errors of different Monte Carlo methods with fitted 
      linear regression lines.
    """
    N_inv_sqrt = N_values**-0.5
    
    methods = [("Monte Carlo", monte_errors), 
                ("Antithetic Variates", anti_errors), 
                ("Moment Matching", moment_errors), 
                ("Halton Sequence", halton_errors)]
    
    plt.figure(figsize=(12, 8))
    
    for label, errors in methods:
        slope, intercept, _, _, _ = linregress(N_inv_sqrt, errors)
        fitted_errors = slope * N_inv_sqrt + intercept
        plt.plot(N_inv_sqrt, errors, 'o', label=f"{label} Errors")
        plt.plot(N_inv_sqrt, fitted_errors, label=f"{label} Fitted Line (y = {slope:.2f}x + {intercept:.2f})")

    plt.xlabel("$N^{-1/2}$")
    plt.ylabel("Error")
    plt.title("Monte Carlo Error vs. $N^{-1/2}$")
    plt.legend()
    plt.grid(True)
    plt.show()

def path_dep_N_analysis():
    """
    Analyzes the behavior of a path-dependent option as the number of simulation 
    paths (N) increases. Plots the option value and errors against N and performs 
    error analysis.

    Returns:
    - Plots showing the path-dependent option value and error as functions of N.
    """
    N_values = np.arange(1000, 20000, 100)
    path_dep_vals = np.zeros(N_values.size)
    path_dep_errors = np.zeros(N_values.size)
    

    for i, N in enumerate(N_values):
        print(i)
        path_dep_vals[i], path_dep_errors[i] = path_dependent_option(N, 30)
        
    mean = np.mean(path_dep_vals)
    print(f"path dependent mean: {mean}")
   
    plt.errorbar(N_values, path_dep_vals, yerr=path_dep_errors, fmt='o', label="Path Dependent Option Value")
    plt.plot(N_values, [mean] * N_values.size, label="Mean Value")
    plt.xlabel("N")
    plt.ylabel("Option Value $C(S_0, t=0)$")
    plt.grid()
    plt.legend()
    plt.show()

    plt.errorbar(N_values, path_dep_errors, fmt='o', label="Path Dependent errors")
    plt.xlabel("N")
    plt.ylabel("error")
    plt.grid()
    plt.legend()
    plt.show()

    error_analysis(N_values, path_dep_errors)

def N_analysis(analytical_value, payoff_function, current_stock_price, expiration_time, interest_rate, mean_function, deviation_function, halton_primes=(None, None)):
    """
    Analyzes the Monte Carlo simulation results as the number of paths (N) increases 
    for different Monte Carlo methods. Compares the simulated option values with 
    the analytical value and performs error analysis.

    Parameters:
    - analytical_value: The known analytical value of the option.
    - payoff_function: The payoff function of the option.
    - current_stock_price: The current stock price (S_0).
    - expiration_time: The expiration time of the option (T).
    - interest_rate: The risk-free interest rate (r).
    - mean_function: The function that computes the mean stock price.
    - deviation_function: The function that computes the standard deviation of the stock price.
    - halton_primes: Tuple of prime numbers for the Halton sequence.

    Returns:
    - Plots comparing the simulated option values with the analytical value and 
      error analysis for different Monte Carlo methods.
    """
    N_values = np.arange(1000, 20000, 100)
    methods = ["normal", "antithetic", "moment_matching", "halton"]
    
    # Initialize a dictionary to store the results for each method
    results = {method: {"values": np.zeros(N_values.size), "errors": np.zeros(N_values.size)} for method in methods}
    
    # Loop through each N value
    for i, N in enumerate(N_values):
        print(f"Processing N = {N} ({i + 1}/{len(N_values)})")
        
        # Loop through each method
        for m in methods:
            # For Halton, we need to pass the halton_primes argument
            if m == "halton":
                results[m]["values"][i], results[m]["errors"][i] = monte_carlo(
                    N, payoff_function, current_stock_price, expiration_time, interest_rate, mean_function, deviation_function, method=m, halton_primes=halton_primes
                )
            else:
                results[m]["values"][i], results[m]["errors"][i] = monte_carlo(
                    N, payoff_function, current_stock_price, expiration_time, interest_rate, mean_function, deviation_function, method=m
                )
    
    print(analytical_value)
    # Plotting
    def plot_results(x, y, yerr=None, xlabel="", ylabel="", label="", fmt='o', title=""):
        plt.errorbar(x, y, yerr=yerr, fmt=fmt, label=label)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid()
        plt.legend()
        plt.title(title)
       

    # Plot 1: Monte Carlo Value vs. Analytical Value
    plot_results(N_values, results["normal"]["values"], results["normal"]["errors"], 
                 xlabel="N", ylabel="Option Value $C(S_0, t=0)$", 
                 label="Monte Carlo Value", title="Monte Carlo vs Analytical Value")
    plt.plot(N_values, [analytical_value] * N_values.size, label="Analytical Value")
    plt.show()

    # Plot 2: Comparison of different methods
    for method in methods:
        plot_results(N_values, results[method]["values"], xlabel="N", ylabel="Option Value $C(S_0, t=0)$", 
                     label=f"{method.capitalize()} Value", fmt="-")
    
    plt.plot(N_values, [analytical_value] * N_values.size, label="Analytical Value")
    plt.xlabel("N")
    plt.ylabel("Option Value $C(S_0, t=0)$")
    plt.grid()
    plt.legend()
    plt.show()

    error_analysis(N_values, results["normal"]["errors"])
    multi_error_analysis(N_values, 
                         results["normal"]["errors"], 
                         results["antithetic"]["errors"], 
                         results["moment_matching"]["errors"], 
                         results["halton"]["errors"])

def time_analysis(num_paths, payoff_function, current_stock_price, expiration_time, interest_rate, mean_function, deviation_function, halton_primes, num_runs=100):
    """
    Analyzes and compares the execution time of different Monte Carlo simulation methods.
    
    Parameters:
    - num_paths: The number of simulation paths (N) to be used in each method.
    - payoff_function: The payoff function of the option.
    - current_stock_price: The current stock price (S_0).
    - expiration_time: The expiration time of the option (T).
    - interest_rate: The risk-free interest rate (r).
    - mean_function: The function that computes the mean.
    - deviation_function: The function that computes the deviation or standard deviation.
    - halton_primes: Tuple of prime numbers for the Halton sequence.
    - num_runs: The number of times each method is executed for timing (default is 100).
    
    Returns:
    - timings: A dictionary containing the timing results for each method.
    """
    methods = ["normal", "antithetic", "moment_matching", "halton"]
    timings = {}

    # Loop through each method
    for m in methods:
        if m == "halton":
            # Time the execution of the Halton sequence method
            time_taken = timeit(
                lambda: monte_carlo(num_paths, payoff_function, current_stock_price, expiration_time, interest_rate, mean_function, deviation_function, method=m, halton_primes=halton_primes),
                number=num_runs
            )
        else:
            # Time the execution of other methods
            time_taken = timeit(
                lambda: monte_carlo(num_paths, payoff_function, current_stock_price, expiration_time, interest_rate, mean_function, deviation_function, method=m),
                number=num_runs
            )
        
        # Store the timing result
        timings[m] = time_taken
        print(f"Time taken for {m} method: {time_taken:.4f} seconds (over {num_runs} runs).")

    return timings

def derivative_analysis(n_sigma):
    derivative_values = np.zeros(n_sigma)
    ds_values = np.zeros(n_sigma)

    for i in range(n_sigma):
        ds = 0.001 + 0.01/n_sigma * i
        ds_values[i] = ds
        derivative_values[i] = path_dependent_derivative(100000, 30, ds)

    ds_slope, ds_intercept, _, _, _ = linregress(ds_values, derivative_values)
    ds_fitted_errors = ds_slope * ds_values + ds_intercept
    print(ds_intercept)

    plt.plot(ds_values, derivative_values, 'o')
    plt.plot(ds_values, ds_fitted_errors, 'r')
    plt.xlabel("dsigma")
    plt.ylabel("Derivative")
    plt.title("Path Dependent Derivative vs. d sigma")
    plt.grid(True)
    plt.show()
