#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cstdlib>
#include <thread>
#include <Windows.h>

std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103,
                           107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241,
                           251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401};

void Choose(const int size, int &first, int &second)
{
    // pick a random element
    first = rand() * size / RAND_MAX;
    // pick a random element from what's left (there is one fewer to choose from)...
    second = rand() * (size - 1) / RAND_MAX;
    // ...and adjust second choice to take into account the first choice
    if (second >= first)
    {
        ++second;
    }
}

double normalDistribution(double x)
{
    // cumulative normal distribution
    if (x < -10.)
        return 0.;
    if (x > 10.)
        return 1.;
    // number of steps
    int N = 2000;
    // range of integration
    double a = 0, b = x;
    // local variables
    double s, h, sum = 0.;
    // inialise the variables
    h = (b - a) / N;
    // add in the first few terms
    sum = sum + exp(-a * a / 2.) + 4. * exp(-(a + h) * (a + h) / 2.);
    // and the last one
    sum = sum + exp(-b * b / 2.);
    // loop over terms 2 up to N-1
    for (int i = 1; i < N / 2; i++)
    {
        s = a + 2 * i * h;
        sum = sum + 2. * exp(-s * s / 2.);
        s = s + h;
        sum = sum + 4. * exp(-s * s / 2.);
    }
    // complete the integral
    sum = 0.5 + h * sum / 3. / sqrt(8. * atan(1.));
    // return result
    return sum;
}

double Halton_Sequence(int i, int b)
{
    // calculates values from the halton sequence
    double f = 1;
    double r = 0;
    while (i > 0)
    {
        f = f / b;
        r = r + f * (i % b);
        i = i / b;
    }
    return r;
}

std::vector<double> box_muller(int a, int b, int N)
{
    // returns normally distributed set of variables from the halton sequence
    std::vector<double> muller_dist(2 * N);
    int k = 1;
    for (int i = 0; i < 2 * N; i = i + 2)
    {
        k++;
        muller_dist[i] = cos(2 * M_PI * Halton_Sequence(k, b)) * sqrt(-2 * log(Halton_Sequence(k, a)));
        muller_dist[i + 1] = sin(2 * M_PI * Halton_Sequence(k, a)) * sqrt(-2 * log(Halton_Sequence(k, b)));
    }
    return muller_dist;
}

std::vector<double> monteCarlo_halton(int M, double S0, double strikePrice, double interestRate, double sigma, double maturity, int N, double alpha, double beta, double theta, double gamma)
{
    // calculates the value returns a vector of values of the option price found from halton sequences monte carlo
    std::vector<double> samples(M);

    for (int i = 0; i < M; i++)
    {
        // declare the distribution
        int a, b;
        Choose(primes.size(), a, b);
        std::vector<double> normal_dist = box_muller(primes[a], primes[b], N);

        // initialise sum
        double sum = 0.;
        double phi;
        double ST;
        for (int i = 0; i < 2 * N; i++)
        {
            phi = normal_dist[i];
            ST = S0 * (alpha * maturity + tan(beta * maturity)) + theta * cosh(2 * beta * maturity - alpha * maturity) + sigma * (1 + beta * maturity) * (1 + beta * maturity) * S0 * exp(2 * gamma - 2 + alpha * maturity) * sqrt(maturity) * phi;
            sum = sum + std::max(strikePrice - ST, 0.0);
        }

        samples[i] = sum / 2 / N * exp(-interestRate * maturity);
    }
    return samples;
}

std::vector<double> monteCarlo_moment(int M, double S0, double strikePrice, double interestRate, double sigma, double maturity, int N, double alpha, double beta, double theta, double gamma)
{
    // calculates the value returns a vector of values of the option price found from moment matching monte carlo
    std::vector<double> samples(M);
    // declare the random number generator
    static std::mt19937 rng;
    // declare the distribution
    std::normal_distribution<> ND(0., 1.);
    ND(rng);

    for (int i = 0; i < M; i++)
    {
        // initialise sum
        double sum = 0.;

        std::vector<double> phi_values(2 * N);
        double vari = 0;

        for (int i = 0; i < N; i = i + 2)
        {
            phi_values[i] = ND(rng);
            phi_values[i + 1] = -phi_values[i];
            vari += 2 * phi_values[i] * phi_values[i];
        }

        vari = vari / N;

        double ST;
        for (int i = 0; i < 2 * N; i++)
        {
            phi_values[i] = phi_values[i] / sqrt(vari);

            ST = S0 * (alpha * maturity + tan(beta * maturity)) + theta * cosh(2 * beta * maturity - alpha * maturity) + sigma * (1 + beta * maturity) * (1 + beta * maturity) * S0 * exp(2 * gamma - 2 + alpha * maturity) * sqrt(maturity) * phi_values[i];
            sum = sum + std::max(strikePrice - ST, 0.0);
        }
        samples[i] = sum / N * exp(-interestRate * maturity);
    }
    return samples;
}

std::vector<double> monteCarlo_antithetic(int M, double S0, double strikePrice, double interestRate, double sigma, double maturity, int N, double alpha, double beta, double theta, double gamma)
{
    // calculates the value returns a vector of values of the option price found from antiothetic variables monte carlo
    std::vector<double> samples(M);
    // declare the random number generator
    static std::mt19937 rng;
    // declare the distribution
    std::normal_distribution<> ND(0., 1.);
    ND(rng);

    for (int i = 0; i < M; i++)
    {
        double phi;
        double STplus;
        double STminus;
        // initialise sum
        double sum = 0.;
        for (int i = 0; i < N; i++)
        {
            phi = ND(rng);
            STplus = S0 * (alpha * maturity + tan(beta * maturity)) + theta * cosh(2 * beta * maturity - alpha * maturity) + sigma * (1 + beta * maturity) * (1 + beta * maturity) * S0 * exp(2 * gamma - 2 + alpha * maturity) * sqrt(maturity) * phi;
            sum = sum + std::max(strikePrice - STplus, 0.0);
            STminus = S0 * (alpha * maturity + tan(beta * maturity)) + theta * cosh(2 * beta * maturity - alpha * maturity) - sigma * (1 + beta * maturity) * (1 + beta * maturity) * S0 * exp(2 * gamma - 2 + alpha * maturity) * sqrt(maturity) * phi;
            sum = sum + std::max(strikePrice - STminus, 0.0);
        }
        samples[i] = sum / 2 / N * exp(-interestRate * maturity);
    }
    return samples;
}

std::vector<double> monteCarlo(int M, double S0, double strikePrice, double interestRate, double sigma, double maturity, int N, double alpha, double beta, double theta, double gamma)
{
    // calculates the value returns a vector of values of the option price found from standard monte carlo
    std::vector<double> samples(M);
    // declare the random number generator
    static std::mt19937 rng;
    // declare the distribution
    std::normal_distribution<> ND(0., 1.);
    ND(rng);

    for (int i = 0; i < M; i++)
    {
        double phi;
        double ST;
        // initialise sum
        double sum = 0.;
        for (int i = 0; i < N; i++)
        {
            phi = ND(rng);
            ST = S0 * (alpha * maturity + tan(beta * maturity)) + theta * cosh(2 * beta * maturity - alpha * maturity) + sigma * (1 + beta * maturity) * (1 + beta * maturity) * S0 * exp(2 * gamma - 2 + alpha * maturity) * sqrt(maturity) * phi;
            sum = sum + std::max(strikePrice - ST, 0.0);
        }
        samples[i] = sum / N * exp(-interestRate * maturity);
    }
    return samples;
}

double P_analytic(double S0, double strikePrice, double interestRate, double sigma, double maturity, double alpha, double beta, double theta, double gamma)
{
    double f = S0 * (alpha * maturity + tan(beta * maturity)) + theta * cosh(2 * beta * maturity - alpha * maturity);
    double v = sigma * (1 + beta * maturity) * (1 + beta * maturity) * S0 * exp(2 * gamma - 2 + alpha * maturity);
    double z = (strikePrice - f) / (v * sqrt(maturity));
    return (strikePrice * normalDistribution(z) + v * sqrt(maturity) * exp(-z * z / 2) / sqrt(2 * M_PI) - f * normalDistribution(z)) * exp(-interestRate * maturity);
}

std::vector<double> data(std::vector<double> samples, double M)
{
    // calculates the value of the mean and variance of normally distributed values
    std::vector<double> params(2);
    double sum = 0.0;
    for (int i = 0; i < M; i++)
    {
        sum += samples[i];
    }
    double mean = sum / M;
    // std::cout << " sample mean = " << mean << std::endl;

    double sumvar = 0.;
    for (int i = 0; i < M; i++)
    {
        sumvar += (samples[i] - mean) * (samples[i] - mean);
    }
    double variance = sumvar / (M - 1);
    // std::cout << " sample variance = " << variance << std::endl;

    double sd = sqrt(variance / M);
    // std::cout << " 95% confident result is in [" << mean - 2. * sd << "," << mean + 2. * sd << "] with " << N * M << " total paths. " << std::endl;

    params[0] = mean;
    params[1] = variance;

    return params;
}

std::vector<double> path_dependent_option(int M, double S0, double strikePrice, double interestRate, double sigma, double maturity, int N, double alpha, double beta, double theta, double gamma, int K)
{
    // calculates the value of the path dependent option
    std::vector<double> samples(M);

    static std::mt19937 rng;
    // declare the distribution
    std::normal_distribution<> ND(0., 1.);
    ND(rng);

    for (int i = 0; i < M; i++)
    {
        std::vector<double> stock_price(K + 1);
        double payoff_sum = 0;
        double sum;
        for (int j = 0; j < N; j++)
        {
            stock_price[0] = S0;
            sum = 0;
            for (int i = 1; i <= K; i++)
            {
                double f = alpha * theta - beta * stock_price[i - 1];
                double v = sigma * pow(fabs(stock_price[i - 1]), gamma);
                double phi = ND(rng);
                stock_price[i] = stock_price[i - 1] + f * maturity / K + v * sqrt(maturity / K) * phi;
                sum += stock_price[i];
            }
            sum = sum / K;
            payoff_sum += std::max(strikePrice - sum, 0.0);
        }
        samples[i] = payoff_sum / N * exp(-interestRate * maturity);
    }
    return samples;
}

std::vector<double> path_dependent_derivative(int M, double S0, double strikePrice, double interestRate, double sigma, double maturity, int N, double alpha, double beta, double theta, double gamma, int K, double d_alpha)
{
    // calculates the derivative of the path dependent option for a certain d_alpha

    std::vector<double> samples_1(M);
    std::vector<double> samples_2(M);

    static std::mt19937 rng;
    // declare the distribution
    std::normal_distribution<> ND(0., 1.);
    ND(rng);

    for (int i = 0; i < M; i++)
    {
        std::vector<double> stock_price_1(K + 1);
        std::vector<double> stock_price_2(K + 1);
        double payoff_sum_1 = 0;
        double payoff_sum_2 = 0;
        double sum_1;
        double sum_2;
        for (int j = 0; j < N; j++)
        {
            stock_price_1[0] = S0;
            stock_price_2[0] = S0;
            sum_1 = 0;
            sum_2 = 0;
            for (int i = 1; i <= K; i++)
            {
                double f_1 = (alpha + d_alpha) * theta - beta * stock_price_1[i - 1];
                double v_1 = sigma * pow(fabs(stock_price_1[i - 1]), gamma);
                double f_2 = (alpha)*theta - beta * stock_price_2[i - 1];
                double v_2 = sigma * pow(fabs(stock_price_2[i - 1]), gamma);
                double phi = ND(rng);
                stock_price_1[i] = stock_price_1[i - 1] + f_1 * maturity / K + v_1 * sqrt(maturity / K) * phi;
                stock_price_2[i] = stock_price_2[i - 1] + f_2 * maturity / K + v_2 * sqrt(maturity / K) * phi;
                sum_1 += stock_price_1[i];
                sum_2 += stock_price_2[i];
            }
            sum_1 = sum_1 / K;
            sum_2 = sum_2 / K;
            payoff_sum_1 += std::max(strikePrice - sum_1, 0.0);
            payoff_sum_2 += std::max(strikePrice - sum_2, 0.0);
        }
        samples_1[i] = payoff_sum_1 / N * exp(-interestRate * maturity);
        samples_2[i] = payoff_sum_2 / N * exp(-interestRate * maturity);
    }

    std::vector<double> derivative_samples(samples_1.size());

    for (int i = 0; i < samples_1.size(); i++)
    {
        derivative_samples[i] = (samples_1[i] - samples_2[i]) / d_alpha;
    }

    return derivative_samples;
}

void function_test()
{
    // calcualtes the mean, standard deviation and confidence interval for every method for a given N
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << "Function Test:" << std::endl;
    // analytic solution
    double P = P_analytic(500.72, 500.0, 0.03, 0.23, 2.0, 0.01, 0.01, 501.0, 1.01);
    std::cout << "analytic value = " << P << std::endl;

    std::vector<double> Monte_normal = monteCarlo(100, 500.72, 500.0, 0.03, 0.23, 2.0, 40000, 0.01, 0.01, 501.0, 1.01);
    std::cout << "normal monteCarlo value = " << data(Monte_normal, Monte_normal.size())[0] << std::endl;
    std::cout << "normal monteCarlo Standard Deviation = " << sqrt(data(Monte_normal, Monte_normal.size())[1] / 100) << std::endl;
    double variance = data(Monte_normal, Monte_normal.size())[1];
    double sd = sqrt(variance / 100);
    std::cout << "95% confident result is in [" << data(Monte_normal, Monte_normal.size())[0] - 2. * sd << "," << data(Monte_normal, Monte_normal.size())[0] + 2. * sd << "] with " << 100000 * 100 << " total paths. " << std::endl;

    std::cout << std::endl;

    std::vector<double> Monte_anti = monteCarlo_antithetic(100, 500.72, 500.0, 0.03, 0.23, 2.0, 40000, 0.01, 0.01, 501.0, 1.01);
    std::cout << "antithetic value = " << data(Monte_anti, Monte_anti.size())[0] << std::endl;
    std::cout << "antithetic Standard Deviation = " << sqrt(data(Monte_anti, Monte_anti.size())[1] / 100) << std::endl;
    variance = data(Monte_anti, Monte_anti.size())[1];
    sd = sqrt(variance / 100);
    std::cout << "95% confident result is in [" << data(Monte_anti, Monte_anti.size())[0] - 2. * sd << "," << data(Monte_anti, Monte_anti.size())[0] + 2. * sd << "] with " << 100000 * 100 << " total paths. " << std::endl;

    std::cout << std::endl;

    std::vector<double> Monte_moment = monteCarlo_moment(100, 500.72, 500.0, 0.03, 0.23, 2.0, 40000, 0.01, 0.01, 501.0, 1.01);
    std::cout << "moment matching value = " << data(Monte_moment, Monte_moment.size())[0] << std::endl;
    std::cout << "moment matching Standard Deviation = " << sqrt(data(Monte_moment, Monte_moment.size())[1] / 100) << std::endl;
    variance = data(Monte_moment, Monte_moment.size())[1];
    sd = sqrt(variance / 100);
    std::cout << "95% confident result is in [" << data(Monte_moment, Monte_moment.size())[0] - 2. * sd << "," << data(Monte_moment, Monte_moment.size())[0] + 2. * sd << "] with " << 100000 * 100 << " total paths. " << std::endl;

    std::cout << std::endl;

    std::vector<double> Monte_halton = monteCarlo_halton(100, 500.72, 500.0, 0.03, 0.23, 2.0, 40000, 0.01, 0.01, 501.0, 1.01);
    std::cout << "halton method value = " << data(Monte_halton, Monte_halton.size())[0] << std::endl;
    std::cout << "halton method Standard Deviation = " << sqrt(data(Monte_halton, Monte_halton.size())[1] / 100) << std::endl;
    variance = data(Monte_halton, Monte_halton.size())[1];
    sd = sqrt(variance / 100);
    std::cout << "95% confident result is in [" << data(Monte_halton, Monte_halton.size())[0] - 2. * sd << "," << data(Monte_halton, Monte_halton.size())[0] + 2. * sd << "] with " << 100000 * 100 << " total paths. " << std::endl;

    std::cout << std::endl;

    // S0, strikePrice, interestRate, sigma, maturity, alpha, beta, theta, gamma, K
    std::vector<double> path_price = path_dependent_option(100, 500.72, 500.0, 0.03, 0.23, 1.0, 40000, 0.01, 0.01, 501.0, 1.01, 40);
    std::cout << "Asian value = " << data(path_price, path_price.size())[0] << std::endl;
    std::cout << "Asian Standard Deviation = " << sqrt(data(path_price, path_price.size())[1] / 100) << std::endl;
    variance = data(path_price, path_price.size())[1];
    sd = sqrt(variance / 100);
    std::cout << "95% confident result is in [" << data(path_price, path_price.size())[0] - 2. * sd << "," << data(path_price, path_price.size())[0] + 2. * sd << "] with " << 100000 * 100 << " total paths. " << std::endl;

    std::cout << std::endl;

    std::cout << "-----------------------------------------------------------" << std::endl;
}

void timing_function()
{
    // calculates the time taken for every method for a values of N
    std::ofstream results_monte;
    results_monte.open("results_timer.csv");
    results_monte << "N"
                  << ", "
                  << "Normal time"
                  << ","
                  << "Antithetic time"
                  << ","
                  << "Moment Matching time"
                  << ","
                  << "Halton time"
                  << ","
                  << "Path Dependent time"
                  << "\n";

    for (int i = 1; i < 5; i++)
    {
        int N = 10000 * i;
        results_monte
            << N
            << ",";

        auto m_StartTime = std::chrono::system_clock::now();
        std::vector<double> Monte_normal = monteCarlo(100, 500.72, 500.0, 0.03, 0.23, 2.0, N, 0.01, 0.01, 501.0, 1.01);
        auto m_endTime = std::chrono::system_clock::now();
        auto m_time = std::chrono::duration_cast<std::chrono::milliseconds>(m_endTime - m_StartTime).count();
        results_monte << m_time / 1000.0
                      << ",";

        auto a_StartTime = std::chrono::system_clock::now();
        std::vector<double> Monte_anti = monteCarlo_antithetic(100, 500.72, 500.0, 0.03, 0.23, 2.0, N, 0.01, 0.01, 501.0, 1.01);
        auto a_endTime = std::chrono::system_clock::now();
        auto a_time = std::chrono::duration_cast<std::chrono::milliseconds>(a_endTime - a_StartTime).count();
        results_monte << a_time / 1000.0
                      << ",";

        auto mm_StartTime = std::chrono::system_clock::now();
        std::vector<double> Monte_moment = monteCarlo_moment(100, 500.72, 500.0, 0.03, 0.23, 2.0, N, 0.01, 0.01, 501.0, 1.01);
        auto mm_endTime = std::chrono::system_clock::now();
        auto mm_time = std::chrono::duration_cast<std::chrono::milliseconds>(mm_endTime - mm_StartTime).count();
        results_monte << mm_time / 1000.0
                      << ",";

        auto h_StartTime = std::chrono::system_clock::now();
        std::vector<double> Monte_halton = monteCarlo_halton(100, 500.72, 500.0, 0.03, 0.23, 2.0, N, 0.01, 0.01, 501.0, 1.01);
        auto h_endTime = std::chrono::system_clock::now();
        auto h_time = std::chrono::duration_cast<std::chrono::milliseconds>(h_endTime - h_StartTime).count();
        results_monte << h_time / 1000.0
                      << ",";

        auto p_StartTime = std::chrono::system_clock::now();
        std::vector<double> path_price = path_dependent_option(100, 500.72, 500.0, 0.03, 0.23, 1.0, N, 0.01, 0.01, 501.0, 1.01, 40);
        auto p_endTime = std::chrono::system_clock::now();
        auto p_time = std::chrono::duration_cast<std::chrono::milliseconds>(p_endTime - p_StartTime).count();
        results_monte << p_time / 1000.0
                      << ","
                      << "\n";
    }
}

void N_range_monte()
{
    // calculates the option value for every method for a range of N values

    // we want to run M calculations
    int M = 100;
    // now store all the results
    std::vector<double> normal_samples(M), anti_samples(M), moment_samples(M), halton_samples(M), path_samples(M);

    std::ofstream results_monte;
    results_monte.open("results_halton.csv");

    results_monte << "N"
                  << ", "
                  << "Normal Mean"
                  << ","
                  << "Normal error"
                  << ","
                  << ","
                  << "Antithetic Mean"
                  << ","
                  << "Antithetic error"
                  << ","
                  << ","
                  << "Moment Matching Mean"
                  << ","
                  << "Moment Matching error"
                  << ","
                  << ","
                  << "Halton Mean"
                  << ","
                  << "Halton error"
                  << ","
                  << ","
                  << "Path Dependent Mean"
                  << ","
                  << "Path Dependent error"
                  << "\n";

    std::cout << "Calculating Monte Carlo simulations for range of N values:" << std::endl;

    for (int k = 1; k < 100; k++)
    {
        int N = 1000 * k;

        normal_samples = monteCarlo(M, 500.72, 500.0, 0.03, 0.23, 2.0, N, 0.01, 0.01, 501.0, 1.01);
        anti_samples = monteCarlo_antithetic(M, 500.72, 500.0, 0.03, 0.23, 2.0, N, 0.01, 0.01, 501.0, 1.01);
        moment_samples = monteCarlo_moment(M, 500.72, 500.0, 0.03, 0.23, 2.0, N, 0.01, 0.01, 501.0, 1.01);
        halton_samples = monteCarlo_halton(M, 500.72, 500.0, 0.03, 0.23, 2.0, N, 0.01, 0.01, 501.0, 1.01);
        path_samples = path_dependent_option(M, 500.72, 500.0, 0.03, 0.23, 1.0, N, 0.01, 0.01, 501.0, 1.01, 40);

        std::vector<double> normal_params = data(normal_samples, M);
        std::vector<double> anti_params = data(anti_samples, M);
        std::vector<double> moment_params = data(moment_samples, M);
        std::vector<double> halton_params = data(halton_samples, M);
        std::vector<double> path_params = data(path_samples, M);

        results_monte << N
                      << ", "
                      << normal_params[0]
                      << ","
                      << sqrt(normal_params[1] / M)
                      << ","
                      << ","
                      << anti_params[0]
                      << ","
                      << sqrt(anti_params[1] / M)
                      << ","
                      << ","
                      << moment_params[0]
                      << ","
                      << sqrt(moment_params[1] / M)
                      << ","
                      << ","
                      << halton_params[0]
                      << ","
                      << sqrt(halton_params[1] / M)
                      << ","
                      << ","
                      << path_params[0]
                      << ","
                      << sqrt(path_params[1] / M)
                      << "\n";

        std::cout << "Current N: " << N << std::endl;
    }

    results_monte.close();

    std::cout << "Results exported to results_monte.csv" << std::endl;
    std::cout << "-----------------------------------------------------------" << std::endl;
}

void derivative_alpha_range()
{
    // calculates the derivative for a range of alpha values
    std::ofstream results_deriv;
    results_deriv.open("results_deriv.csv");

    int N = 10000;

    results_deriv << "N"
                  << ", "
                  << "d_alpha"
                  << ","
                  << "Derivative"
                  << "\n";

    for (int i = 1; i < 20; i++)
    {
        double d_alpha = i * 0.0001;
        std::cout << d_alpha << std::endl;
        std::vector<double> deriv = path_dependent_derivative(100, 500.72, 500.0, 0.03, 0.23, 1.0, 30000, 0.01, 0.01, 501.0, 1.01, 40, d_alpha);
        results_deriv << N
                      << ", "
                      << d_alpha
                      << ","
                      << data(deriv, deriv.size())[0]
                      << "\n";
    }
    results_deriv.close();
}

int main()
{
    // include functions for tests needing to be taken.
    function_test();
    timing_function();
    N_range_monte();
    derivative_alpha_range();
    return 0;
}