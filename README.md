# ODE_approximation

$\textit{Included in this repository is:}$
1. A program to simulate particle motion using Runge-Kutta
2. Interfaces for numerically approximating ODEs
3. Approximations of the Runge-Kutta coefficients for any order using a genetic algorithm

$\textbf{Dependencies}$

Your going to need OpenCV and OpenMP for this one.

$\textbf{Particle simulation:}$
This was done in order for me to qualitatively compare the accuracy of different Runge-Kutta orders.

In the "particle_simulation.cpp" file, the user may change the force field by modifying the f_x and f_y functions.

The program will ask for:
1. Time step size (h)
2. Number of time steps
3. The order of Runge-Kutta approximation, and the text file containing the coefficients
4. The number of particles, and for each particle it's initial condition.

It will then show the path of the particle(s) at normal speed.

$\textbf{General numerical approximations simulations:}$

The "Approximations.h" file contains an interface called Difeq_Sim that, given a defined system and IVP, can calculate a table of values of the system over a certain time interval. 

Later versions of this project may include different approximations such as Newton's forward interpolation or Taylor series.

$\textbf{Approximations of the Runge-Kutta coefficients with a GA:}$

Runge-Kutta coefficients are hard to obtain for any general order. The "ApproximateCoefficients.cpp" file frames this as an optimization problem.

A genetic algorithm is an optimization technique that simulates the process of natural selection by maintaining a population where each member is ranked on how well it fits the ideal condition.

In this case, our chromosones are defined by the Runge-Kutta coefficients.

Goal function:
Minimize the error between the true solution and predicted values of the solution governed by a chromosone's RK coefficients. 

Results:
Sadly, I was only able to meet the local truncation error requirements for RK-2. The coefficients for higher orders would not give a terribly far off approximation, however.

