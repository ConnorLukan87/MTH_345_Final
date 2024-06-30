//
// Created by Connor Lukan on 4/20/24.
//

#ifndef RUNGEKUTTAMTH345_TYPES_H
#define RUNGEKUTTAMTH345_TYPES_H

#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <math.h>


typedef long long int ll;

using Function = std::function<double(std::vector<double>)>;

namespace maths {
    const int factorial[7] = {1, 1, 2, 6, 24, 120, 720};

    int triangle(int i)
    {
        return i*(i+1)/2;
    }

    std::pair<int, int> transform(double& x, double& y, double& m, int scale)
    {
        int col = ((double)scale)*(.5*x/m + .5) - 1.0;
        int row = ((double)scale)*(1.0-(.5*y/m + .5));
        return std::make_pair(row, col);
    }

    double theta(double x1, double y1, double x2, double y2)
    {
        return std::atan((y2-y1)/(x2-x1));
    }

    double distance(double x1, double y1, double x2, double y2)
    {
        return std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    }

};

class coefficients {
public:
    std::vector<double> a, v, b;

    coefficients(int N)
    {
        a = std::vector<double>(N+1);
        b = std::vector<double>((N-1)*N/2);
        v = std::vector<double>(N+1);
    }

    coefficients() {}

    static coefficients load_coefficients(int N, std::string filename)
    {
        coefficients temp_coeffs = coefficients(N);
        std::ifstream coefficients_file(filename);
        if (!coefficients_file.is_open())
        {
            throw std::invalid_argument("Yo dawg. File couldn't open for some reason.");
        }

        std::string temp = "";
        for (int v_j = 1; v_j <= N; v_j++)
        {
            coefficients_file >> temp;
            temp_coeffs.v[v_j] = std::stod(temp);
        }

        // read b
        int sz = temp_coeffs.b.size();
        for (int b_k = 0 ; b_k < sz; b_k++)
        {
            coefficients_file >> temp;
            temp_coeffs.b[b_k] = std::stod(temp);
        }

        // read a

        sz = temp_coeffs.a.size();
        for (int a_l = 2; a_l <= N; a_l++)
        {
            coefficients_file >> temp;
            temp_coeffs.a[a_l] = std::stod(temp);
        }

        coefficients_file.close();
        return temp_coeffs;
    }
};

struct particle
{
    double m;

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> v_x;
    std::vector<double> v_y;
    std::vector<double> t;

    particle(int steps_total)
    {
        x = std::vector<double>(steps_total+1);
        y = std::vector<double>(steps_total+1);
        v_x = std::vector<double>(steps_total+1);
        v_y = std::vector<double>(steps_total+1);
        t = std::vector<double>(steps_total+1);
    }

    particle()
    {;}
};

typedef struct {
    std::vector<double> t;
    std::vector<std::vector<double>> x;
} Path;



#endif //RUNGEKUTTAMTH345_TYPES_H
