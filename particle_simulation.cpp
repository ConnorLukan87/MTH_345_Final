//
// Created by Connor Lukan on 4/20/24.
//

#include <iostream>
#include <filesystem>
#include <opencv2/opencv.hpp>
#include "types.h"

namespace fs = std::filesystem;

const int OUTPUT_RESOLUTION = 512;
std::string PATH = "/home/computer/ODE_approximation/";


double f_x(double x, double y);
double f_y(double x, double y);
double f_inter(int state_i, particle& p1, particle& p2, int direction); // interparticle force
double s_0_prime(std::vector<double>& params); // velocity
double s_1_prime(std::vector<double>& params, int state_i, int particle_id, std::vector<particle> all_particles, int direction); // acceleration. computed w Newtons second.
void update(int state_i, int particle_id, std::vector<particle>& all_particles, double& h, coefficients& coeffs); // step for 1 particle

int main()
{
    std::cout << "Step size: ";
    double h;
    std::cin >> h;

    int k_0;
    std::cout << "How many steps: ";
    std::cin >> k_0;

    std::cout << "Which approximation would you like to use: \n0) Runge-Kutta\n";

    int choice;
    std::cin >> choice;

    std::cout << "How many particles: ";
    int m;
    std::cin >> m;

    std::vector<particle> particles(m, particle(k_0));

    double max_bounds = 0.0;

    for (int p_i = 0; p_i < m ; p_i++)
    {
        std::cout << "PARTICLE " << p_i << ":\n";
        std::cout << "Mass: ";
        std::cin >> particles[p_i].m;
        std::cout << "Initial position (x, y): ";
        std::cin >> particles[p_i].x[0] >> particles[p_i].y[0];
        std::cout << "Initial velocity (x, y): ";
        std::cin >> particles[p_i].v_x[0] >> particles[p_i].v_y[0];
        particles[p_i].t[0] = 0.0;

        std::cout << std::endl;
    }

    if (choice == 0)
    {
        // rk
        std::cout << "What order runge-kutta: ";
        int N;
        std::cin >> N;
        std::cout << "Coefficients file name: ";
        std::string coeffs_file;
        std::cin >> coeffs_file;
        coefficients coeffs = coefficients::load_coefficients(N, PATH + coeffs_file);

        // set the min max bounds
        for (int p = 0 ; p < m ; p++)
        {
            max_bounds = std::max(max_bounds, std::abs(particles[p].x[0]));
            max_bounds = std::max(max_bounds, std::abs(particles[p].y[0]));
        }

        for (int state_i = 0 ; state_i < k_0; state_i++)
        {
            for (int particle_j = 0 ; particle_j < m; particle_j++)
            {
                update(state_i, particle_j, particles, h, coeffs);
                double x_pos = particles[particle_j].x[state_i+1];
                double y_pos = particles[particle_j].y[state_i+1];

                if (std::abs(x_pos) > max_bounds)
                {
                    max_bounds = std::abs(x_pos);
                }

                if (std::abs(y_pos) > max_bounds)
                {
                    max_bounds = std::abs(y_pos);
                }
            }
        }

    }
    else
    {
        std::cout << "Invalid choice. Exiting.\n";
    }

    cv::Mat last(OUTPUT_RESOLUTION, OUTPUT_RESOLUTION, CV_8UC1);

    cv::namedWindow("Simulated particle motion", cv::WINDOW_NORMAL);
    cv::resizeWindow("Simulated particle motion", cv::Size(OUTPUT_RESOLUTION, OUTPUT_RESOLUTION));



    for (int t_i = 0; t_i < k_0; t_i++)
    {

        cv::Mat to_show(last);

        for (int particle_j = 0; particle_j < particles.size(); particle_j++)
        {

            std::pair<int,int> coords = maths::transform(particles[particle_j].x[t_i], particles[particle_j].y[t_i], max_bounds, OUTPUT_RESOLUTION);

            if (coords.second < 0)
            {
                coords.second = 0;
            }

            if (coords.first < 0)
            {
                coords.first = 0;
            }

            if (coords.first >= OUTPUT_RESOLUTION)
            {
                coords.first = OUTPUT_RESOLUTION - 1;
            }

            if (coords.second >= OUTPUT_RESOLUTION)
            {
                coords.second = OUTPUT_RESOLUTION - 1;
            }

            to_show.at<uchar>(coords.first, coords.second) = 255;
        }

        cv::imshow("Simulated particle motion", to_show);

        int q = cv::waitKey(1000.0*h); // 1000ms/second * h seconds

        if (q == 'q') // exit
        {
            cv::destroyAllWindows();
            break;
        }
        last = to_show;
    }

    cv::waitKey(3000); // leave it up for 3 seconds
    cv::destroyAllWindows();

    int amount = 0;


    // write it to a file
    std::string path = PATH + "RK_simulation_outs/";

    for (const auto& file : fs::directory_iterator(path))
    {
        amount++;
    }

    std::ofstream out_file(path + "out_" + std::to_string(amount) + ".txt", std::ios::out);

    if (!out_file.is_open())
    {
        std::cout << "YO: Output file could not open." << std::endl;
        return -1;
    }

    for (int j = 0 ; j < m; j++)
    {
        // write mass and particle id
        out_file << j << "\n" << particles[j].m << std::endl;
        for (int i= 0 ; i <= k_0; i++) // for every step
        {
            out_file << particles[j].t[i] << "\t" << particles[j].x[i] << "\t" << particles[j].y[i] << "\t" << particles[j].v_x[i] << "\t" << particles[j].v_y[i] << std::endl;
        }
        out_file << std::endl;
    }

    out_file.close();
    return 0;

}

double f_x(double x, double y)
{
    return 0.0;
}

double f_y(double x, double y)
{
    return 0.0;//-1.0*x + (1.0-(x*x) - (y*y))*y;
}


double f_inter(int state_i, particle& p1, particle& p2, int direction)
{
    double force = p1.m*p2.m;

    if (p1.x[state_i] < p2.x[state_i] && p1.y[state_i] < p2.y[state_i])
    {
        double angle = maths::theta(p1.x[state_i], p1.y[state_i], p2.x[state_i], p2.y[state_i]);
        double dist = maths::distance(p1.x[state_i], p1.y[state_i], p2.x[state_i], p2.y[state_i]);
        force = force/(dist*dist);
        return force*(direction == 0 ? std::cos(angle) : std::sin(angle));
    }
    else if (p1.x[state_i] >= p2.x[state_i] && p1.y[state_i] >= p2.y[state_i])
    {
        return -f_inter(state_i, p2, p1, direction);
    }
    else if (p1.x[state_i] >= p2.x[state_i] && p1.y[state_i] <= p2.y[state_i])
    {
        double angle = maths::theta(p1.x[state_i], p1.y[state_i], p2.x[state_i], p2.y[state_i]) + M_PI;
        double dist = maths::distance(p1.x[state_i], p1.y[state_i], p2.x[state_i], p2.y[state_i]);
        force = force/(dist*dist);
        return force*(direction == 0 ? std::cos(angle) : std::sin(angle));
    }
    else if (p1.x[state_i] < p2.x[state_i] && p1.y[state_i] > p2.y[state_i])
    {
        return -f_inter(state_i, p2, p1, direction);
    }


   /* double angle = maths::theta(p1.x[state_i], p1.y[state_i], p2.x[state_i], p2.y[state_i]);
    double dist = maths::distance(p1.x[state_i], p1.y[state_i], p2.x[state_i], p2.y[state_i]);
    force = force/(dist*dist);
    return force*(direction == 0 ? std::cos(angle) : std::sin(angle));*/
}


double s_0_prime(std::vector<double>& params)
{
    return params[2];
}

double s_1_prime(std::vector<double>& params, int state_i, int particle_id, std::vector<particle> all_particles, int direction)
{
    double rv = direction == 0 ? f_x(all_particles[particle_id].x[state_i], all_particles[particle_id].y[state_i]) : f_y(all_particles[particle_id].x[state_i], all_particles[particle_id].y[state_i]);

    for (int i = 0 ; i < all_particles.size(); i++)
    {
        if (i != particle_id)
        {
            rv += f_inter(state_i, all_particles[particle_id], all_particles[i], direction);
        }
    }

    return rv/all_particles[particle_id].m;
}


void update(int state_i, int particle_id, std::vector<particle>& all_particles, double& h, coefficients& coeffs)
{

    for (int direction = 0 ; direction < 2; direction++)
    {
        std::vector<std::vector<double>> f(2, std::vector<double>(5));

        std::vector<double> g(2);

        std::vector<double> temp(3);
        temp[0] = all_particles[0].t[state_i];
        temp[1] = direction == 0 ? all_particles[particle_id].x[state_i] : all_particles[particle_id].y[state_i];
        temp[2] = direction == 0 ? all_particles[particle_id].v_x[state_i] : all_particles[particle_id].v_y[state_i];

        f[0][1] = s_0_prime(temp);
        f[1][1] = s_1_prime(temp, state_i, particle_id, all_particles, direction);

        for (int n = 2; n <= 4; n++)
        {
            for (int k_j = 0 ; k_j < 2; k_j++)
            {
                std::vector<double> parameters(3);
                parameters[0] = all_particles[0].t[state_i] + h*coeffs.a[n];

                for (int p = 1; p <= 2; p++)
                {
                    for (int l = 1; l < n; l++)
                    {
                        parameters[p] += h*coeffs.b[maths::triangle(n-2) + l-1]*f[p-1][l];
                    }
                }

                f[k_j][n] = (k_j == 0 ? s_0_prime(parameters) : s_1_prime(parameters, state_i, particle_id, all_particles, direction));
            }
        }

        if (direction == 0)
        {
            all_particles[particle_id].x[state_i+1] = 0;
            all_particles[particle_id].v_x[state_i+1] = 0;
            for (int v_i = 1 ; v_i <= 4; v_i++)
            {
                all_particles[particle_id].x[state_i+1] += f[0][v_i]*coeffs.v[v_i];
                all_particles[particle_id].v_x[state_i+1] += f[1][v_i]*coeffs.v[v_i];
            }
            all_particles[particle_id].x[state_i+1] *= h;
            all_particles[particle_id].v_x[state_i+1] *= h;

            all_particles[particle_id].x[state_i+1] += all_particles[particle_id].x[state_i];
            all_particles[particle_id].v_x[state_i+1] += all_particles[particle_id].v_x[state_i];

        }
        else
        {
            all_particles[particle_id].y[state_i+1] = 0;
            all_particles[particle_id].v_y[state_i+1] = 0;
            for (int v_i = 1 ; v_i <= 4; v_i++)
            {
                all_particles[particle_id].y[state_i+1] += f[0][v_i]*coeffs.v[v_i];
                all_particles[particle_id].v_y[state_i+1] += f[1][v_i]*coeffs.v[v_i];
            }
            all_particles[particle_id].y[state_i+1] *= h;
            all_particles[particle_id].v_y[state_i+1] *= h;

            all_particles[particle_id].y[state_i+1] += all_particles[particle_id].y[state_i];
            all_particles[particle_id].v_y[state_i+1] += all_particles[particle_id].v_y[state_i];
        }
    }

}
