//
// Created by Connor Lukan on 4/20/24.
//

#ifndef RUNGEKUTTAMTH345_APPROXIMATIONS_H
#define RUNGEKUTTAMTH345_APPROXIMATIONS_H
#include "types.h"

class Difeq_Simulation {
protected:
    int K; // k many equations
    std::vector<Function> x_prime;

    double t_0;
    std::vector<double> x_0;
public:
    Difeq_Simulation(int K, std::vector<Function> x_prime, double t_0, std::vector<double> x_0)
    {
        this->t_0 = t_0;
        this->x_0 = x_0;
        this->x_prime = x_prime;
        this->K = K;
    }

    Difeq_Simulation() = delete;

    virtual Path solve_ivp(double t_f);
};

class RK_Simulation : public Difeq_Simulation {
private:
    int N;
    coefficients coeffs;
    double h;
    bool N_set;

public:
    RK_Simulation(int k, std::vector<Function>& x_prime, int t_0, std::vector<double> x_0, int N, double h) : Difeq_Simulation(k, x_prime, t_0, x_0)
    {
        this->N = N;
        N_set = true;
        this->h = h;

        this->setCoefficients("/home/computer/CS/RungeKuttaMTH345/default_RK_coeffs/RK_" + std::to_string(N) + ".txt");
    }


    void setCoefficients(std::string filename)
    {
        this->coeffs = coefficients::load_coefficients(N, filename);
    }

    Path solve_ivp(double t_final)
    {

        if (t_final < t_0)
        {
            throw std::invalid_argument("Final time less than initial time and time machine is unavailable.");
        }
        ll total_steps = (t_final-t_0)/h;

        Path path;
        // path.t = std::vector<double>(total_steps+1);
        path.x = std::vector<std::vector<double>>();
        path.x.push_back(x_0);

        std::vector<std::vector<double>> f(K, std::vector<double>(N+1));

        double t_i = t_0;
        std::vector<double> current_x = x_0;

        std::vector<double> temp;
        for (ll i = 0; i < total_steps; i++) {
            std::vector<double> g(K); // g[i] will contain the running sum of x[i](t_0 + h)
            // set these to x(t_i))

            // f[][1] vector will be x_prime(t_i)

            temp = std::vector<double>(K + 1);
            temp[0] = t_i;
            for (int k = 1; k <= K; k++)
            {
                temp[k] = current_x[k-1];
            }

            // prepare matrix for dynamic programming
            for (int j = 0 ; j < K; j++)
            {
                f[j][1] = x_prime[j](temp);
            }

            for (int n = 2 ; n <= N; n++)
            {
                for (int k_j = 0 ; k_j < K; k_j++)
                {
                    // calculate this f_n_j

                    std::vector<double> parameters(K+1);
                    parameters[0] = t_i + h*this->coeffs.a[n];

                    // for each of the other parameters, do the summation of the beta_n values
                    for (int p = 1 ; p <= K; p++)
                    {
                        // use the p^th row of the f vector and get all f_p_i's for i in range(1, n)
                        parameters[p] = current_x[p-1];
                        for (int l = 1 ; l < n; l++)
                        {
                            parameters[p] += h*this->coeffs.b[maths::triangle(n-2) + l-1]*f[p-1][l]; // ahh yeah this is confusing sorry
                        }
                    }

                    f[k_j][n] = x_prime[k_j](parameters); // its the evaluation of x_prime[k_j]
                }
            }

            // now, just evaluate! everthing is in the f[][] matrix
            for (int equation = 0 ; equation < K; equation++)
            {
                g[equation] = 0;
                for (int v_i = 1; v_i <=N; v_i++)
                {
                    g[equation] += f[equation][v_i]*this->coeffs.v[v_i];
                }
                g[equation] *= h;
                g[equation] += current_x[equation];
            }

            t_i += h;
            current_x = g;
            path.t.push_back(t_i);
            path.x.push_back(current_x);
        }

        return path;
    }

};



#endif //RUNGEKUTTAMTH345_APPROXIMATIONS_H
