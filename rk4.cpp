#include <iostream>
#include <cmath>
#include <valarray>
#include "constants.hpp"

using namespace std;


// Runge-Kutta 4th order method
valarray<double> RungeKutta4(double t0,
                             valarray<double> y0,
                             double h,
                             double tf) {
    double t = t0;
    valarray<double> y = y0;

    while (t < tf) {
        valarray<double> k1 = h * rhs(t, y);
        valarray<double> k2 = h * rhs(t + h / 2.0, y + k1 / 2.0);
        valarray<double> k3 = h * rhs(t + h / 2.0, y + k2 / 2.0);
        valarray<double> k4 = h * rhs(t + h, y + k3);

        // Update y using the weighted average of slopes
        y += (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;

        // Move to the next step
        t += h;
    }

    return y;
}


