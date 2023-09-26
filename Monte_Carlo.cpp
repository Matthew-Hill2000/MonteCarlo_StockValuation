#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <algorithm>
#include <fstream>
#include <vector>
#include <functional>
#include <cstdlib>
#include <ctime>
#include <chrono>
using namespace std;

double cumNormDis(double x)
{
    return 0.5 * erfc(-x / sqrt(2.));
}

///////////////////////////////////////// Monte-Carlo Functions /////////////////////////////////////

double monteCarlo(double S0, double strikePrice, double interestRate, double maturity, int N,
                  double alpha, double beta, double theta, double volatility, double gamma)
{
    // declare the random number generator
    static mt19937 rng;
    // declare the distribution
    normal_distribution<> ND(0., 1.);
    // initialise sum
    double sum = 0.;
    for (int i = 0; i < N; i++)
    {
        double phi = ND(rng);
        double ST = S0 * (alpha * maturity + tan(beta * maturity)) + theta * cos(alpha * maturity + beta * maturity) + volatility * (1 + alpha * maturity) * 0.5 * pow(S0 + theta, gamma) * pow(maturity, 0.5) * phi;
        sum = sum + max(strikePrice - ST, 0.);
    }
    return sum / N * exp(-interestRate * maturity);
}

double monteCarlo_antithetic(double S0, double strikePrice, double interestRate, double maturity, int N,
                             double alpha, double beta, double theta, double volatility, double gamma)
{
    // declare the random number generator
    static mt19937 rng;
    // declare the distribution
    normal_distribution<> ND(0., 1.);
    // initialise sum
    double sum = 0.;
    for (int i = 0; i < N; i++)
    {
        double phi = ND(rng);
        double ST_plus = S0 * (alpha * maturity + tan(beta * maturity)) + theta * cos(alpha * maturity + beta * maturity) + volatility * (1 + alpha * maturity) * 0.5 * pow(S0 + theta, gamma) * pow(maturity, 0.5) * phi;
        sum = sum + max(strikePrice - ST_plus, 0.);
        double ST_minus = S0 * (alpha * maturity + tan(beta * maturity)) + theta * cos(alpha * maturity + beta * maturity) + volatility * (1 + alpha * maturity) * 0.5 * pow(S0 + theta, gamma) * pow(maturity, 0.5) * -phi;
        sum = sum + max(strikePrice - ST_minus, 0.);
    }
    return sum / 2 / N * exp(-interestRate * maturity);
}

double monteCarlo_with_vector(double S0, double strikePrice, double interestRate, double maturity,
                              double alpha, double beta, double theta, double volatility, double gamma, std::vector<double> phi_vector)
{
    // initialise sum
    double sum = 0.;
    int N = phi_vector.size();
    for (int i = 0; i < N; i++)
    {
        double phi = phi_vector[i];
        double ST = S0 * (alpha * maturity + tan(beta * maturity)) + theta * cos(alpha * maturity + beta * maturity) + volatility * (1 + alpha * maturity) * 0.5 * pow(S0 + theta, gamma) * pow(maturity, 0.5) * phi;
        sum = sum + max(strikePrice - ST, 0.);
    }
    return sum / N * exp(-interestRate * maturity);
}

////////////////////////////////////// Vector function //////////////////////////////////////////////////

std::vector<double> phi_vector_func(int p1, int p2, int N)
{
    std::vector<double> phi_vector(2 * N);
    static mt19937 rng;
    normal_distribution<> ND(0., 1.);
    double phi;
    double sum{0};

    for (int i = 0; i < N; i++)
    {
        phi = ND(rng);
        phi_vector[2 * i] = phi;
        phi_vector[2 * i + 1] = -phi;
        sum += 2 * phi * phi;
    }

    sum /= (2 * N - 1);
    double std{sqrt(sum)};
    for (int i = 0; i < 2 * N; i++)
    {
        phi_vector[i] /= std;
    }
    return phi_vector;
}
////////////////////////////////// Halton Formulas //////////////////////////////////////////////

double Halton_Seq(int index, int base)
{
    double f = 1, r = 0;
    while (index > 0)
    {
        f = f / base;
        r = r + f * (index % base);
        index = index / base;
    }
    return r;
}

std::vector<std::vector<double>> Halton_vector_func(int a, int b, int N)
{
    std::vector<std::vector<double>> Halton_vector(N, vector<double>(2));
    for (int i = 0; i < N; i++)
    {
        Halton_vector[i][0] = Halton_Seq(i + 1, a);
        Halton_vector[i][1] = Halton_Seq(i + 1, b);
    }
    return Halton_vector;
}

std::vector<double> box_muller_func(std::vector<std::vector<double>> Halton_vector)
{
    std::vector<double> box_muller_vec(2 * Halton_vector.size());
    for (int i = 0; i < Halton_vector.size(); i++)
    {
        double x1 = Halton_vector[i][0];
        double x2 = Halton_vector[i][1];
        // double r = sqrt(pow(x1, 2) + pow(x2, 2));
        box_muller_vec[2 * i] = cos(2 * atan(1) * 4 * x2) * sqrt(-2 * log(x1));     //= x1*sqrt((-2 * log(r)) / r);
        box_muller_vec[2 * i + 1] = sin(2 * atan(1) * 4 * x1) * sqrt(-2 * log(x2)); //= x2*sqrt((-2 * log(r)) / r);
    }
    return box_muller_vec;
}

std::vector<double> Halton_phi_vector_func(int p1, int p2, int N)
{
    return box_muller_func(Halton_vector_func(p1, p2, N));
}

///////////////////////////Generating random primes/////////////////////////////////////

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

//////////////////////// Analytic Formulas //////////////////////////////////////////

double f_func(double S0, double strikePrice, double interestRate, double maturity, int N,
              double alpha, double beta, double theta, double volatility, double gamma)
{
    return S0 * (alpha * maturity + tan(beta * maturity)) + theta * cos(alpha * maturity + beta * maturity);
}

double v_func(double S0, double strikePrice, double interestRate, double maturity, int N,
              double alpha, double beta, double theta, double volatility, double gamma)
{
    return volatility * (1 + alpha * maturity) * 0.5 * pow(S0 + theta, gamma);
}

double z_func(double S0, double strikePrice, double interestRate, double maturity, int N,
              double alpha, double beta, double theta, double volatility, double gamma)
{
    return (strikePrice - f_func(S0, strikePrice, interestRate, maturity, N,
                                 alpha, beta, theta, volatility, gamma)) /
           (v_func(S0, strikePrice, interestRate, maturity, N,
                   alpha, beta, theta, volatility, gamma) *
            pow(maturity, 0.5));
}

double put_func(double S0, double strikePrice, double interestRate, double maturity, int N,
                double alpha, double beta, double theta, double volatility, double gamma)
{
    return ((strikePrice * cumNormDis(z_func(S0, strikePrice, interestRate, maturity, N,
                                             alpha, beta, theta, volatility, gamma))) +
            v_func(S0, strikePrice, interestRate, maturity, N,
                   alpha, beta, theta, volatility, gamma) *
                pow(maturity / (2. * atan(1) * 4), 0.5) * exp(-0.5 * pow(z_func(S0, strikePrice, interestRate, maturity, N, alpha, beta, theta, volatility, gamma), 2)) -
            f_func(S0, strikePrice, interestRate, maturity, N, alpha, beta, theta, volatility, gamma) *
                cumNormDis(z_func(S0, strikePrice, interestRate, maturity, N, alpha, beta, theta, volatility, gamma))) *
           exp(-interestRate * maturity);
}

//////////////////////////////////////Output to file functions////////////////////////////////////////////////////////////

void output_to_file(double S0, double strikePrice, double interestRate, double maturity, int N, double alpha, double beta, double theta,
                    double volatility, double gamma, function<double(double S0, double strikePrice, double interestRate, double maturity, int N, double alpha, double beta, double theta, double volatility, double gamma)> func, std::string file_name, int step_size, int max_N)
{
    ofstream myfile;
    myfile.open(file_name);
    for (int i = 1; i <= max_N / step_size; i++)
    {
        myfile << step_size * i << "," << func(S0, strikePrice, interestRate, maturity, step_size * i, alpha, beta, theta, volatility, gamma) << "\n";
    }
    myfile.close();
}

void output_to_file_vector(double S0, double strikePrice, double interestRate, double maturity, int N, double alpha, double beta, double theta,
                           double volatility, double gamma, int p1, int p2, function<double(double S0, double strikePrice, double interestRate, double maturity, double alpha, double beta, double theta, double volatility, double gamma, std::vector<double> phi_vector)> mc_func, function<std::vector<double>(int p1, int p2, int N)> vector_func, std::string file_name, int step_size, int max_N)
{
    ofstream myfile;
    myfile.open(file_name);
    for (int i = 1; i <= max_N / step_size; i++)
    {
        std::vector<double> phi_vector = vector_func(p1, p2, step_size * i);
        myfile << step_size * i << "," << mc_func(S0, strikePrice, interestRate, maturity, alpha, beta, theta, volatility, gamma, phi_vector) << "\n";
    }
    myfile.close();
}

///////////////////////////////////////////////Using Monte Carlo//////////////////////////////////////////////////////////////

void average_errors(double S0, double strikePrice, double interestRate, double maturity, int N, double alpha, double beta, double theta,
                    double volatility, double gamma, int step_size, int max_N, int repeats, std::string file_name, double solution)
{

    ofstream myfile;
    myfile.open(file_name);
    myfile << "N"
           << ",";
    for (int i = 1; i <= max_N / step_size; i++)
    {
        myfile << step_size * i << ",";
    }
    myfile << "\n";

    for (int i = 1; i <= repeats; i++)
    {
        myfile << "error " << i << ",";
        for (int j = 1; j <= max_N / step_size; j++)
        {
            double error = fabs(monteCarlo(S0, strikePrice, interestRate, maturity, step_size * j,
                                           alpha, beta, theta, volatility, gamma) -
                                solution);
            myfile << error << ",";
        }
        myfile << "\n";
    }
    myfile.close();
}

std::vector<double> average_montecarlo(double S0, double strikePrice, double interestRate, double maturity, int N, double alpha, double beta, double theta,
                                       double volatility, double gamma, int repeats, function<double(double S0, double strikePrice, double interestRate, double maturity, int N, double alpha, double beta, double theta, double volatility, double gamma)> func)
{
    vector<double> samples(repeats);

    for (int i = 0; i < repeats; i++)
    {
        samples[i] = func(S0, strikePrice, interestRate, maturity, N,
                          alpha, beta, theta, volatility, gamma);
    }
    // estimate the mean from the sample
    double sum = 0.;
    for (int i = 0; i < repeats; i++)
    {
        sum += samples[i];
    }
    double mean = sum / repeats;
    // estimate the variance from the sample
    double sumvar = 0.;
    for (int i = 0; i < repeats; i++)
    {
        sumvar += (samples[i] - mean) * (samples[i] - mean);
    }
    double variance = sumvar / (repeats - 1);

    // get the standard deviation of the sample mean
    double sd = sqrt(variance / repeats);

    std::vector<double> mean_and_error = {mean, sd};
    return mean_and_error;
}

std::vector<double> average_montecarlo_vector(double S0, double strikePrice, double interestRate, double maturity, int N,
                                              double alpha, double beta, double theta, double volatility, double gamma, int repeats, int p1, int p2, function<std::vector<double>(int p1, int p2, int N)> vector_func, function<double(double S0, double strikePrice, double interestRate, double maturity, double alpha, double beta, double theta, double volatility, double gamma, std::vector<double> vector)> mc_func)
{
    std::vector<double> samples(repeats);

    for (int i = 0; i < repeats; i++)
    {
        std::vector<double> phi_vector = vector_func(p1, p2, N);
        samples[i] = mc_func(S0, strikePrice, interestRate, maturity,
                             alpha, beta, theta, volatility, gamma, phi_vector);
    }
    // estimate the mean from the sample
    double sum = 0.;
    for (int i = 0; i < repeats; i++)
    {
        sum += samples[i];
    }
    double mean = sum / repeats;
    // estimate the variance from the sample
    double sumvar = 0.;
    for (int i = 0; i < repeats; i++)
    {
        sumvar += (samples[i] - mean) * (samples[i] - mean);
    }
    double variance = sumvar / (repeats - 1);

    // get the standard deviation of the sample mean
    double sd = sqrt(variance / repeats);

    std::vector<double> mean_and_error = {mean, sd};

    return mean_and_error;
}

std::vector<double> average_montecarlo_Halton_vector(double S0, double strikePrice, double interestRate, double maturity, int N,
                                                     double alpha, double beta, double theta, double volatility, double gamma, int repeats, function<std::vector<double>(int p1, int p2, int N)> Halton_phi_vector_func, function<double(double S0, double strikePrice, double interestRate, double maturity, double alpha, double beta, double theta, double volatility, double gamma, std::vector<double> vector)> mc_func)
{
    std::vector<double> samples(repeats);
    int a;
    int b;

    for (int i = 0; i < repeats; i++)
    {
        Choose(primes.size(), a, b);
        int p1 = primes[a];
        int p2 = primes[b];
        std::vector<double> Halton_phi_vector = Halton_phi_vector_func(p1, p2, N);
        samples[i] = mc_func(S0, strikePrice, interestRate, maturity,
                             alpha, beta, theta, volatility, gamma, Halton_phi_vector);
    }
    // estimate the mean from the sample
    double sum = 0.;
    for (int i = 0; i < repeats; i++)
    {
        sum += samples[i];
    }
    double mean = sum / repeats;
    // estimate the variance from the sample
    double sumvar = 0.;
    for (int i = 0; i < repeats; i++)
    {
        sumvar += (samples[i] - mean) * (samples[i] - mean);
    }
    double variance = sumvar / (repeats - 1);

    // get the standard deviation of the sample mean
    double sd = sqrt(variance / repeats);

    std::vector<double> mean_and_error = {mean, sd};

    return mean_and_error;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main()
{
    /*double valueOld = 1., diffOld = 1.; int n = 100;
    for (int i = 1; i <= 14; i++)
    {
        n *= 2;
        // calculate value with n
        double value = monteCarlo(70.0231, 70., 0.01, 1., n, 0.02, 0.02, 70., 0.19, 1.03);
        // and difference from last time
        double diff = fabs(value - valueOld);
        // output stage, steps, value,
        cout << i << " " << n << " " << value << " " << diff << endl;

        // store old values
        valueOld = value;
        diffOld = diff;
    }*/

    // output_to_file(70.0231, 70., 0.01, 1., 10000, 0.02, 0.02, 70., 0.19, 1.03, &monteCarlo_antithetic, "2.1.3_data_antithetic.csv", 500, 50000);
    // output_to_file_vector(70.0231, 70., 0.01, 1., 10000, 0.02, 0.02, 70., 0.19, 1.03, 2, 3, &monteCarlo_with_vector, &phi_vector_func,
    //"2.1.3_data_momentmatch.csv", 500, 50000);
    // output_to_file_vector(70.0231, 70., 0.01, 1., 10000, 0.02, 0.02, 70., 0.19, 1.03, 241, 2, &monteCarlo_with_vector, &Halton_phi_vector_func,
    //"2.1.3_data_Halton.csv", 500, 50000);

    /*ofstream myfile_2;
myfile_2.open("2.1.3_data_antithetic.csv");
for (int i = 1; i <= 1000; i++){
myfile_2 << 100*i << "," << monteCarlo_antithetic(70.0231, 70., 0.01, 1., 100*i, 0.02, 0.02, 70., 0.19, 1.03) << "\n";}
myfile_2.close();*/

    /*std::vector<double> Bm = Halton_phi_vector_func(2, 3, 5000);
    std::vector<double> phis{ phi_vector_func(10000) };
    std::cout << monteCarlo_with_vector(70.0231, 70., 0.01, 1., 0.02, 0.02, 70., 0.19, 1.03, Bm) << std::endl;*/

    double anal_soln = put_func(70.0231, 70., 0.01, 1., 10000, 0.02, 0.02, 70., 0.19, 1.03);

    // average_errors(70.0231, 70., 0.01, 1., 10000, 0.02, 0.02, 70., 0.19, 1.03, 500, 50000, 100, "2.1.2_convergence.csv", anal_soln);

    // std::vector<double> a = average_montecarlo(70.0231, 70., 0.01, 1., 10000, 0.02, 0.02, 70., 0.19, 1.03, 1000, &monteCarlo);
    // std::vector<double> answer = average_montecarlo_Halton_vector(70.0231, 70., 0.01, 1., 5000, 0.02, 0.02, 70., 0.19, 1.03, 1000,
    //&Halton_phi_vector_func, &monteCarlo_with_vector);

    /*ofstream myfile_3;
    myfile_3.open("2.1.3_accuracy_and_times.csv");
    myfile_3 << "Method" << "," << "mean_value" << "," << "mean_error" << "," << "time_taken" << "\n";

    auto antithetic_start = std::chrono::steady_clock::now();
    std::vector<double> antithetic_me = average_montecarlo(70.0231, 70., 0.01, 1., 5000, 0.02, 0.02, 70., 0.19, 1.03, 1000, &monteCarlo_antithetic);
    auto antithetic_finish = std::chrono::steady_clock::now();
    auto antithetic_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(antithetic_finish - antithetic_start);
    myfile_3 << "Antithetic" << "," << antithetic_me[0] << "," << antithetic_me[1] << "," << antithetic_elapsed.count() << "\n";

    auto momentmatch_start = std::chrono::steady_clock::now();
    std::vector<double> momentmatch_me = average_montecarlo_vector(70.0231, 70., 0.01, 1., 5000, 0.02, 0.02, 70., 0.19, 1.03, 1000, &phi_vector_func,
        &monteCarlo_with_vector);
    auto momentmatch_finish = std::chrono::steady_clock::now();
    auto momentmatch_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(momentmatch_finish - momentmatch_start);
    myfile_3 << "Moment Match" << "," << momentmatch_me[0] << "," << momentmatch_me[1] << "," << momentmatch_elapsed.count() << "\n";

    auto Halton_start = std::chrono::steady_clock::now();
    std::vector<double> Halton_me = average_montecarlo_Halton_vector(70.0231, 70., 0.01, 1., 5000, 0.02, 0.02, 70., 0.19, 1.03, 1000, &Halton_phi_vector_func,
        &monteCarlo_with_vector);
    auto Halton_finish = std::chrono::steady_clock::now();
    auto Halton_elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(Halton_finish - Halton_start);
    myfile_3 << "Halton" << "," << Halton_me[0] << "," << Halton_me[1] << "," << Halton_elapsed.count() << "\n";

    myfile_3.close();*/

    return 0;
}