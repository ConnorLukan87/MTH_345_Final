//
// Created by Connor Lukan on 4/20/24.
//

#include "types.h"
#include <math.h>

/*
 * This file verifies that the RK-N coefficients satisfy the local truncation error requirement
 */

std::string PATH = "/home/computer/CS/RungeKuttaMTH345/";

void load_population(std::vector<coefficients>& population, std::string& path, int& N);
bool test_coefficients(coefficients& coeffs, double t_0, double x_0, double h, double steps, int& N);
double x_eval(coefficients& chromo, double x_i, double t_i, double& h, int& N);
double x(double t, double t_i, double x_i);
double x_nplus1(double t, double x);
double g(double t_i, double x_i);

int main()
{
    std::cout << "Which order RK: ";
    int N;
    std::cin >> N;

    std::cout << "Which output file: ";

    int out_i;
    std::cin >> out_i;

    std::string path = PATH + "RK_" + std::to_string(N) + "_Coeffs/out_" + std::to_string(out_i) + "/";

    std::cout << "How many members of the population would you like to load: ";
    int sz;
    std::cin >> sz;

    std::vector<coefficients> population(sz, coefficients(N));

    load_population(population, path, N);

    bool ans = test_coefficients(population[0], 0.0, 1.0, 0.01, 1000, N);
    if (ans)
    {
        std::cout << "Hooray! The test passed.";
    }
    else
    {
        std::cout << "Oh no! Coefficients are not sufficient!" << std::endl;
    }
}

void load_population(std::vector<coefficients>& population, std::string& path, int& N)
{
    std::ifstream coefficients_file(path + "coeffs.txt");

    if (!coefficients_file.is_open())
    {
        std::cout << "Coefficinets file could not open!" << std::endl;
        return;
    }

    std::cout << "Population size: " << population.size() << std::endl;
    for (int member_i = 0 ; member_i < population.size(); member_i++)
    {
        // read v
        std::string temp = "";
        for (int v_j = 1; v_j <= N; v_j++)
        {
            coefficients_file >> temp;
            population[member_i].v[v_j] = std::stod(temp);
        }

        // read b
        int sz = population[member_i].b.size();
        for (int b_k = 0 ; b_k < sz; b_k++)
        {
            coefficients_file >> temp;
            population[member_i].b[b_k] = std::stod(temp);
        }

        // read a

        sz = population[member_i].a.size();
        for (int a_l = 2; a_l <= N; a_l++)
        {
            coefficients_file >> temp;
            population[member_i].a[a_l] = std::stod(temp);
        }

        coefficients_file >> temp >> temp >> temp;
        std::cout << "Temp:" << std::endl;
    }
    coefficients_file.close();
}

double g(double t_i, double x_i)
{
    return std::exp(t_i) + x_i*x_i;
}

double x_nplus1(double t, double x)
{
    double x_prime = g(t, x);
    double x_2prime = std::exp(t) + 2*x*x_prime;
    double x_3prime = std::exp(t) + 2*(std::pow(g(t,x), 2.0) + x*x_2prime);
    double x_4prime = std::exp(t) + 2*(3*x_prime*x_2prime + x*x_3prime);
    return std::exp(t) + 2*(3*x_2prime*x_2prime + 2*x_prime*x_3prime + x*x_4prime);
}

double x(double t, double t_i, double x_i)
{
    double x = x_i;
    double x_prime = g(t_i, x_i);
    x += (t-t_i)*x_prime;
    double x_2prime = std::exp(t_i) + 2*x_i*x_prime;
    x += std::pow(t-t_i, 2.0)*(x_2prime)/((double)maths::factorial[2]);
    double x_3prime = std::exp(t_i) + 2*(std::pow(g(t_i,x_i), 2.0) + x_i*x_2prime);
    x += std::pow(t-t_i, 3.0)*(x_3prime)/((double)maths::factorial[3]);
    double x_4prime = std::exp(t_i) + 2*(3*x_prime*x_2prime + x_i*x_3prime);
    x += std::pow(t-t_i, 4.0)*(x_4prime)/((double)maths::factorial[4]);
    double x_5prime = std::exp(t_i) + 2*(3*x_2prime*x_2prime + 2*x_prime*x_3prime + x_i*x_4prime);
    x += std::pow(t-t_i, 5.0)*x_5prime/((double)maths::factorial[5]);
    double x_6prime = std::exp(t_i) + 2*(8*x_2prime*x_3prime + 3*x_prime*x_4prime + x_5prime*x_i);
    x += std::pow(t-t_i, 6.0)*x_6prime/((double)maths::factorial[6]);
    return x;
}

double x_eval(coefficients& chromo, double x_i, double t_i, double& h, int& N)
{
    std::vector<double> f_i(N+1);
    f_i[1] = g(t_i, x_i);
    double sum = x_i + chromo.v[1]*h*f_i[1];
    for (int i = 2 ; i < N+1; i++)
    {
        double x_component = x_i;
        for (int b_i_j = 0 ; b_i_j < i-1; b_i_j++)
        {
            x_component += h*chromo.b[(((i-1)*(i-2))/(2)) + b_i_j]*f_i[b_i_j + 1];
        }
        f_i[i] = g(t_i + chromo.a[i]*h, x_component);
        sum += chromo.v[i]*h*f_i[i];
    }
    return sum;
}

double max_xnplus1(double& t_1, double x_1, double t_2)
{
    int intervals = 24;

    double step_size = (t_2-t_1)/((double)intervals);

    double max_val = x_nplus1(t_1, x_1);
    double last_x = x_1;
    for (int i = 1 ; i <= intervals; i++)
    {
        double this_x = x(t_1 + step_size*i, t_1 + step_size*(i-1), last_x);
        double ynplus1 = x_nplus1(t_1 + step_size*i, this_x);
        if (ynplus1 > max_val)
        {
            max_val = ynplus1;
        }
        last_x = this_x;
    }

    return max_val;
}

bool test_coefficients(coefficients& coeffs, double t_0, double x_0, double h, double steps, int& N)
{

    double t_i = t_0;
    double x_i = x_0;
    for (int step_i = 0 ; step_i < steps; step_i++)
    {
        double max = max_xnplus1(t_i, t_i + h, N);
        double actual = x(t_i + h, t_i, x_i);
        double upper_error = (max*(std::pow(h, (double)(N+1))))/((maths::factorial[N+1]));
        double predicted = x_eval(coeffs, x_i, t_i, h, N);
        std::cout << "Max on [" << t_i << ", " << t_i + h << "] = " << max << "\n";
        std::cout << "Error allowed: " << upper_error << std::endl;
        std::cout << "True error: " << std::abs(predicted - actual) << std::endl;
        std::cout << "Predicted: " << predicted << " \tTrue: " << actual << "\n\n";
        if (upper_error < std::abs(predicted - actual)) // );
        {
            return false;
        }
        t_i += h;

    }
    return true;
}
