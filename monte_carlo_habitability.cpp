// monte_carlo_habitability.cpp
// Monte Carlo estimate of the probability that each TRAPPIST-1 planet lies
// in the habitable zone, treating L*, T*, M_p and r as independent Gaussian
// random variables sampled from the means and 1-sigma uncertainties in
// the TRAPPIST-1 data in constants.hpp. The d_HZ_inner / d_HZ_outer / is_habitable functions are
// reproduced here so this file builds as a standalone program.

#include <iostream>
#include <random>
#include <cmath>
#include "constants.hpp"

// Inner habitable-zone boundary in AU (Kopparapu et al. habitable-zone fit).
double d_HZ_inner(double starLum, double starTemp, double planetMass) {
    double dT = starTemp - 5780.0;
    double mEarthUnits = planetMass / M_EARTH;

    double poly = d_inner_coeffs[1] * dT
                + d_inner_coeffs[2] * dT * dT
                + d_inner_coeffs[3] * dT * dT * dT
                + d_inner_coeffs[4] * dT * dT * dT * dT;
    double Seff = d_inner_coeffs[0] + poly * pow(mEarthUnits, 0.1);

    return sqrt(starLum / Seff);
}

// Outer habitable-zone boundary in AU.
double d_HZ_outer(double starLum, double starTemp, double planetMass) {
    double dT = starTemp - 5780.0;
    double mEarthUnits = planetMass / M_EARTH;

    double poly = d_outer_coeffs[1] * dT
                + d_outer_coeffs[2] * dT * dT
                + d_outer_coeffs[3] * dT * dT * dT
                + d_outer_coeffs[4] * dT * dT * dT * dT;
    double Seff = d_outer_coeffs[0] + poly * pow(mEarthUnits, 0.1);

    return sqrt(starLum / Seff);
}

// Returns true if the orbital radius (AU) is between the inner and outer HZ.
bool is_habitable(double starLum, double starTemp,
                  double planetMass, double orbitalRadius) {
    double dInner = d_HZ_inner(starLum, starTemp, planetMass);
    double dOuter = d_HZ_outer(starLum, starTemp, planetMass);
    return (orbitalRadius >= dInner && orbitalRadius <= dOuter);
}

// Helper for runs ns Monte Carlo samples for one planet and returns
// the fraction of samples that fall in the habitable zone.
double estimate_probability(double meanMass, double sigmaMass,
                            double meanRadius, double sigmaRadius,
                            long ns) {
    long pinhz = 0;
    std::mt19937 gen(std::random_device{}());

    // Star properties shared across all planets
    std::normal_distribution<double> disL(L_1a, dL_1a);
    std::normal_distribution<double> disT(T_1a, dT_1a);
    // Planet-specific properties
    std::normal_distribution<double> dism(meanMass, sigmaMass);
    std::normal_distribution<double> disr(meanRadius, sigmaRadius);

    for (long i = 0; i < ns; ++i) {
        if (is_habitable(disL(gen), disT(gen), dism(gen), disr(gen))) pinhz++;
    }
    return static_cast<double>(pinhz) / ns;
}

// ---------------- Single-planet estimate ----------------
int main() {
    long ns = 10000000;     // number of Monte Carlo samples
    long pinhz = 0;

    // Gaussian distributions for each measurement, using the median and
    // 1-sigma uncertainty from constants.hpp.
    std::mt19937 gen(std::random_device{}());
    std::normal_distribution<double> dism(M_1e, dM_1e);
    std::normal_distribution<double> disr(OR_1e_AU, dOR_1e_AU);
    std::normal_distribution<double> disL(L_1a, dL_1a);
    std::normal_distribution<double> disT(T_1a, dT_1a);

    // Sample all four parameters per iteration; count how often the planet
    // lands inside the habitable zone.
    for (long i = 0; i < ns; ++i) {
        double M = dism(gen);
        double r = disr(gen);
        double L = disL(gen);
        double T = disT(gen);
        if (is_habitable(L, T, M, r)) pinhz++;
    }

    double prob = static_cast<double>(pinhz) / ns;
    std::cout << "Estimated using " << ns << " points: " << prob;

    // ---------------- All planets ----------------
    // Repeat the Monte Carlo for every TRAPPIST-1 planet.
    std::cout << "\nHabitability probabilities for all planets\n";
    std::cout << "TRAPPIST-1b: " << estimate_probability(M_1b, dM_1b, OR_1b_AU, dOR_1b_AU, ns) << std::endl;
    std::cout << "TRAPPIST-1c: " << estimate_probability(M_1c, dM_1c, OR_1c_AU, dOR_1c_AU, ns) << std::endl;
    std::cout << "TRAPPIST-1d: " << estimate_probability(M_1d, dM_1d, OR_1d_AU, dOR_1d_AU, ns) << std::endl;
    std::cout << "TRAPPIST-1e: " << estimate_probability(M_1e, dM_1e, OR_1e_AU, dOR_1e_AU, ns) << std::endl;
    std::cout << "TRAPPIST-1f: " << estimate_probability(M_1f, dM_1f, OR_1f_AU, dOR_1f_AU, ns) << std::endl;
    std::cout << "TRAPPIST-1g: " << estimate_probability(M_1g, dM_1g, OR_1g_AU, dOR_1g_AU, ns) << std::endl;
    std::cout << "TRAPPIST-1h: " << estimate_probability(M_1h, dM_1h, OR_1h_AU, dOR_1h_AU, ns) << std::endl;
    return 0;
}