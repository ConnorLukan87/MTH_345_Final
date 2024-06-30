//
// Created by Connor Lukan on 4/20/24.
//

#include <random>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include <filesystem>
#include "types.h"

// these need to be dialed somehow
#define ITERATIONS 1000
#define PERCENT_TO_DIE_OFF 0.25
#define STEP_SIZE 0.0001
#define MUTATION_RATE 0.6
#define MUTATION_SIZE 0.3
#define K_ROYAL 10
#define POPULATION_MULTIPLE 1000
#define STEPS_PER_EVOLUTION 1000

namespace fs = std::filesystem;

const int N = 4; // Runge-Kutta order

std::string PATH = "/home/computer/CS/RungeKuttaMTH345/";

coefficients randomize(coefficients c, std::mt19937& gen);
double max_xnplus1(double a, double b, double x);
void sort(std::vector<std::pair<coefficients, double>>& population);
void score_on(std::vector<std::pair<coefficients, double>>& population, std::vector<double>& y, std::vector<double>& E, double& t_0, double& h);
void mutate(std::vector<std::pair<coefficients, double>>& population, std::mt19937& gen);
double eval(coefficients& chromo, double y_i, double t_i, double& h);
double x(double t, double t_i, double x_i);
double x_nplus1(double t, double x);
double g(double t, double y);

int main()
{
    int n = STEPS_PER_EVOLUTION;
    double h = STEP_SIZE;

    // input initial condiiton. This must be on the predefined Function
    std::cout << "Initial time: ";
    double t_0;
    std::cin >> t_0;

    double y_0;
    std::cout << "Initial y: ";
    std::cin >> y_0;

    std::random_device rd;
    std::mt19937 gen(rd());

    int pop_size = (2*N + (N*(N-1)/2))*POPULATION_MULTIPLE;

    std::vector<std::pair<coefficients, double>> population;

    for (int i = 0 ; i < pop_size; i++)
    {
        population.push_back(std::make_pair(randomize(coefficients(N), gen), 0.0));
    }

    std::vector<double> all_y_calculated, all_E_calculated;
    std::uniform_int_distribution<> av(0, N-1); // split between drawn and drawn +1
    std::uniform_int_distribution<> b(0,-2 + (N-1)*N/2);

    double t_initial = t_0;
    double last = 0.0;
    double last_y = y_0;

    for (int k = 1 ; k <= ITERATIONS; k++)
    {
        auto begin_time = std::chrono::high_resolution_clock::now();
        std::vector<double> y_calculated(n+1), E(n);
        y_calculated[0] = last_y;
        // fill y vector

#pragma omp parallel for
        for (int i = 1; i < n+1; i++)
        {
            y_calculated[i] = x(t_0 + h*i, t_0 + (i-1)*h, y_calculated[i-1]);
        }

#pragma omp parallel for
        for (int i = 0 ; i < n; i++)
        {
            // calculate max error for interval i, i+1
            E[i] = max_xnplus1(t_0 + i*h, t_0 + (i+1)*h, y_calculated[i])*std::pow(h, N+1)/((double)maths::factorial[N+1]);
        }


        score_on(population, y_calculated, E, t_0, h); // score based on new iteration

        // sort
        sort(population);

        int to_die = PERCENT_TO_DIE_OFF*pop_size; // get rid of the worst "to_die" many
        for (int i = 0 ; i < to_die; i++)
        {
            population.pop_back();
        }


        // reproduce, then run the entire simulation on the offspring including this interval
        std::vector<std::pair<coefficients, double>> to_add(to_die, std::pair<coefficients, double>(coefficients(N), 0.0));
        int sz = population.size();

        std::uniform_int_distribution<int> random_other_chromosone(K_ROYAL, sz-1);

        std::cout << "h";
        for (int i = 0; i < to_die; i++)
        {
            coefficients first = population[i%K_ROYAL].first, second = population[random_other_chromosone(gen)].first;

            int v_split = av(gen);
            int j;
            for (j = 0 ; j < v_split; j++)
            {
                to_add[i].first.v[j] = first.v[j];
            }

            for (; j < N+1; j++)
            {
                to_add[i].first.v[j] = second.v[j];
            }


            j = 0;
            int a_split = av(gen);
            for (; j < a_split; j++)
            {
                to_add[i].first.a[j] = first.a[j];
            }

            for (; j < N+1; j++)
            {
                to_add[i].first.a[j] = second.a[j];
            }

            if (N > 2)
            {
                int b_split = b(gen);
                j = 0;
                for (; j < b_split; j++)
                {
                    to_add[i].first.b[j] = first.b[j];
                }

                for (;j < (N-1)*N/2; j++)
                {
                    to_add[i].first.b[j] = second.b[j];
                }
            }
            else
            {
                to_add[i].first.b[0] = i % 2 == 0 ? first.b[0] : second.b[0];//+ second.b[0])/2.0; // maybe do a lil something here
            }

        }
        t_0 = t_0 + h*n;

        for (int i = 1 ; i < N+1; i++)
        {
            all_y_calculated.push_back(y_calculated[i]); // add these y_calculated values to the long y vector, starting at y_calcluated[1]
            all_E_calculated.push_back(E[i-1]);  // add these error values to the long error vector
        }

        score_on(to_add, all_y_calculated, all_E_calculated, t_initial, h);
        // evaluate the offspring on all of these, score them and stuff


        // add ths offspring to the population
        for (auto& member : to_add)
        {
            population.push_back(member);
        }

        sort(population);

        mutate(population, gen);


        auto time_taken = std::chrono::high_resolution_clock::now() - begin_time;

        std::cout << "Population first score: " << population[0].second << std::endl << "\t\tDifference: " << population[0].second - last << std::endl;
        std::cout << "Time taken: " << (double)(std::chrono::duration_cast<std::chrono::milliseconds>(time_taken).count())/1000.0 << " seconds" << std::endl;
        last = population[0].second;
        if (last > std::pow(10.0, 49.0))
        {
            // STOP! Numbers are too big );
            break;
        }
        last_y = all_y_calculated[all_y_calculated.size()-1];
    }

    score_on(population, all_y_calculated, all_E_calculated, t_initial, h);
    sort(population);
    // print out coefficients for the best 3, then compare to the desired result (a.k.a governing equations)
    int amount = 0;


    for (const auto& file : fs::directory_iterator(PATH+ "RK_" + std::to_string(N) + "_Coeffs/"))
    {
        amount++;
    }

    std::string path = PATH + "RK_" + std::to_string(N) + "_Coeffs/out_" + std::to_string(amount) + "/";

    system(("mkdir " + path).c_str());

    std::ofstream params_file(path + "algorithm_params.txt");

    params_file << ITERATIONS << std::endl << STEPS_PER_EVOLUTION << std::endl << STEP_SIZE << std::endl << MUTATION_RATE << std::endl << MUTATION_SIZE << std::endl << K_ROYAL << std::endl << POPULATION_MULTIPLE << std::endl;

    params_file.close();

    std::ofstream coefficients_file(path + "coeffs.txt");

    for (int i = 0 ; i < 30; i ++)
    {
        std::cout << i << " score: " << population[i].second << std::endl;
        std::cout << "v: ";
        for (int j = 1; j <= N;j++)
        {
            std::cout << population[i].first.v[j] << " ";
            coefficients_file << population[i].first.v[j] << " ";
        }
        std::cout << std::endl << "b: \n";
        coefficients_file << std::endl;

        for (int j=0; j < population[i].first.b.size(); j++)
        {
            std::cout << population[i].first.b[j] << " ";
            coefficients_file << population[i].first.b[j] << " ";
        }
        std::cout << std::endl;
        coefficients_file << std::endl;

        std::cout << "a: ";
        for (int j = 2; j < N+1; j++)
        {
            std::cout << population[i].first.a[j] << " ";
            coefficients_file << population[i].first.a[j] << " ";
        }
        std::cout << "\n\n";

        coefficients_file << "\n\n\n";
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

//double x(double t, double )

double eval(coefficients& chromo, double y_i, double t_i, double& h)
{
    std::vector<double> f_i(N+1);
    f_i[1] = g(t_i, y_i);
    double sum = y_i + chromo.v[1]*h*f_i[1];
    for (int i = 2 ; i < N+1; i++)
    {
        double y_component = y_i;
        for (int b_i_j = 0 ; b_i_j < i-1; b_i_j++)
        {
            y_component += h*chromo.b[maths::triangle(i-2) + b_i_j]*f_i[b_i_j + 1];
        }
        f_i[i] = g(t_i + chromo.a[i]*h, y_component);
        sum += chromo.v[i]*h*f_i[i];
    }
    return sum;
}

void mutate(std::vector<std::pair<coefficients, double>>& population, std::mt19937& gen)
{
    int amount = MUTATION_RATE*population.size();
    std::uniform_int_distribution<> rand_pop_size(0, population.size()-1), rand_av(0, N), rand_b(0, -1 + ((N-1)*N/2));
    std::uniform_real_distribution<> coef(-1.0,1.0);
    for (int i = 0 ; i < amount; i++)
    {
        int member_i = rand_pop_size(gen);
        int v_mutation = rand_av(gen);
        int a_mutation = rand_av(gen);
        int b_mutation = rand_b(gen);

        population[member_i].first.v[v_mutation] += coef(gen)*MUTATION_SIZE;
        population[member_i].first.a[a_mutation] += coef(gen)*MUTATION_SIZE;
        population[member_i].first.b[b_mutation] += coef(gen)*MUTATION_SIZE;
    }
}

void score_on(std::vector<std::pair<coefficients, double>>& population, std::vector<double>& y, std::vector<double>& E, double& t_0, double& h) {
    int amount_of_tests = E.size();
    int threads_wanted = 28;
    int amount_per_thread = amount_of_tests / threads_wanted;
    int remainder = amount_of_tests - (threads_wanted * amount_per_thread);

    omp_lock_t pop_second_lock;
    omp_init_lock(&pop_second_lock);
#pragma omp parallel num_threads(threads_wanted)
    {
        int thread_idx = omp_get_thread_num();
        int beginE_idx = thread_idx * amount_per_thread;
        int end_idx = beginE_idx + amount_per_thread;
        double t = t_0 + beginE_idx * h;
        std::vector<double> changes(population.size());
        for (int i = beginE_idx; i < end_idx; i++) {
            for (int j = 0; j < population.size(); j++) {
                double out = eval(population[j].first, y[i], t, h);
                if (out < y[i + 1]) {
                    if (y[i + 1] - out > E[i]) {
                        changes[j] += y[i + 1] - out - E[i];
                    }
                } else {
                    if (out - y[i + 1] > E[i]) {
                        changes[j] += out - (y[i + 1] + E[i]);
                    }
                }
            }
            t += h;
        }
        omp_set_lock(&pop_second_lock);
        for (int j = 0; j < population.size(); j++) {
            population[j].second += changes[j];
        }
        omp_unset_lock(&pop_second_lock);

    }
    omp_destroy_lock(&pop_second_lock);
    for (int i = amount_of_tests - remainder; i < amount_of_tests; i++) { // rest of the interval
#pragma omp parallel for
        for (int j = 0; j < population.size(); j++) {
            double out = eval(population[j].first, y[i], t_0 + i * h, h);
            if (out < y[i + 1]) {
                if (y[i + 1] - out > E[i]) {
                    population[j].second += y[i + 1] - out - E[i];
                }
            } else {
                if (out - y[i + 1] > E[i]) {
                    population[j].second += out - (y[i + 1] + E[i]);
                }
            }
        }
    }

    // add error due to the sum of v_i's = 1 lemma
    // boost it with a constant
#pragma omp parallel for
    for (int i = 0; i < population.size(); i++) {
        double sum_of_v = 0.0;
        for (int j = 1; j <= N; j++) {
            sum_of_v += population[i].first.v[j];
        }
        double sum_of_b;
        double a_val;

        double eqn_err = std::abs(sum_of_v - 1.0);

        for (int j = 1; j < N; j++) {
            a_val = population[i].first.a[j + 1];
            sum_of_b = 0.0;
            for (int l = maths::triangle(j - 1); l < maths::triangle(j); l++) {
                sum_of_b += population[i].first.b[l];
            }

            eqn_err += std::abs(sum_of_b - a_val) - 3.0 * (a_val < 0.0 ? a_val : 0.0);
        }

        population[i].second += 5.0 * eqn_err;
    }

}

void sort(std::vector<std::pair<coefficients, double>>& population)
{
    std::sort(population.begin(), population.end(), [&](std::pair<coefficients, double> c1, std::pair<coefficients, double> c2) {
        return c1.second < c2.second;
    });

}

coefficients randomize(coefficients c, std::mt19937& gen)
{

    std::uniform_real_distribution<double> uniform(-1.0,1.0);

    int i = 0;
    for (;i < N+1; i++)
    {
        c.a[i] = uniform(gen);
        c.v[i] = uniform(gen);
        c.b[i] = uniform(gen);
    }

    for (;i < (N-1)*N/2; i++)
    {
        c.b[i] = uniform(gen);
    }

    return c;
}

double max_xnplus1(double a, double b, double x_0)
{
    // compute the max of the N+1th derivative of y on [a,b]
    int R = 10;
    // maybe do R many points on the interval
    double amount = (b-a)/(double)R;

    double max = std::abs(x_nplus1(a, x_0));
    double t_i = a;
    double x_i =x_0;

    for (int i = 1; i <= R; i++)
    {
        x_i = x(t_i + amount, t_i, x_i); //
        t_i += amount;

        double d = std::abs(x_nplus1(t_i, x_i));
        if (d > max)
        {
            max = d;
        }
    }
    return max;
}
