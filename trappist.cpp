#include <iostream>
#include <iomanip>
#include <fstream>
#include <valarray>
#include <vector>
#include <list>
#include <cmath>
#include <cstdlib>
#include "constants.hpp"
#include "Planet.hpp"
#include "System.hpp"

using namespace std;

// Global system so rhs_NBody can access the masses
System trappist_system(
    "TRAPPIST-1", "TRAPPIST-1a",
    0.0, 0.0, 0.0, 0.0,
    valarray<double>{0.0, 0.0},
    valarray<double>{0.0, 0.0},
    {}
);

// ============================================================
// General N-body RHS
// State vector layout:
// [x0, y0, vx0, vy0, x1, y1, vx1, vy1, ...]
// body 0 = star, bodies 1..N = planets
// ============================================================
valarray<double> rhs_NBody(double t, const valarray<double> yvec) {
    (void)t;

    int nbodies = static_cast<int>(trappist_system.planetList.size()) + 1;
    valarray<double> dydt(4 * nbodies);

    // Collect masses into an indexable vector
    vector<double> masses(nbodies);
    masses[0] = trappist_system.star_mass;

    int p = 1;
    for (const Planet& planet : trappist_system.planetList) {
        masses[p++] = planet.mass;
    }

    for (int i = 0; i < nbodies; i++) {
        double xi = yvec[4*i];
        double yi = yvec[4*i + 1];

        // Position derivatives = velocities
        dydt[4*i]     = yvec[4*i + 2];
        dydt[4*i + 1] = yvec[4*i + 3];

        // Velocity derivatives = gravitational acceleration
        double ax = 0.0;
        double ay = 0.0;

        for (int j = 0; j < nbodies; j++) {
            if (i == j) continue;

            double dx = xi - yvec[4*j];
            double dy = yi - yvec[4*j + 1];
            double r2 = dx*dx + dy*dy;
            double r3 = pow(r2, 1.5);

            ax -= G_SI * masses[j] * dx / r3;
            ay -= G_SI * masses[j] * dy / r3;
        }

        dydt[4*i + 2] = ax;
        dydt[4*i + 3] = ay;
    }

    return dydt;
}

// Wrapper used by rk4.cpp
valarray<double> rhs(double t, const valarray<double> yvec) {
    return rhs_NBody(t, yvec);
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <h> <N>" << endl;
        return 1;
    }

    double h = atof(argv[1]);
    int N = atoi(argv[2]);
    double t = 0.0;

    // ========================================================
    // Build TRAPPIST-1 system using constants from constants.hpp
    // ========================================================
    valarray<double> starX = {0.0, 0.0};
    valarray<double> starV = {0.0, 0.0};

    const int nplanets = 7;

    const string planet_names[nplanets] = {
        "TRAPPIST-1b", "TRAPPIST-1c", "TRAPPIST-1d",
        "TRAPPIST-1e", "TRAPPIST-1f", "TRAPPIST-1g", "TRAPPIST-1h"
    };

    const double radii_AU[nplanets] = {
        OR_1b_AU, OR_1c_AU, OR_1d_AU,
        OR_1e_AU, OR_1f_AU, OR_1g_AU, OR_1h_AU
    };

    const double planet_masses[nplanets] = {
        M_1b, M_1c, M_1d, M_1e, M_1f, M_1g, M_1h
    };

    list<Planet> planets;

    // Place planets at different angular positions around the star
    for (int k = 0; k < nplanets; k++) {
        double r = radii_AU[k] * AU_IN_M;
        double v = sqrt(G_SI * M_1a / r);
        double theta = 2.0 * acos(-1.0) * k / nplanets;

        valarray<double> X = {r * cos(theta), r * sin(theta)};
        valarray<double> V = {-v * sin(theta), v * cos(theta)};

        planets.push_back(Planet(planet_names[k], planet_masses[k], 0.0, X, V));
    }

    trappist_system = System(
        "TRAPPIST-1", "TRAPPIST-1a",
        M_1a, R_1a, T_1a, L_1a,
        starX, starV, planets
    );

    // ========================================================
    // Build state vector from the system
    // ========================================================
    int nbodies = static_cast<int>(trappist_system.planetList.size()) + 1;
    valarray<double> y(4 * nbodies);

    // Star = body 0
    y[0] = trappist_system.star_X[0];
    y[1] = trappist_system.star_X[1];
    y[2] = trappist_system.star_V[0];
    y[3] = trappist_system.star_V[1];

    int i = 1;
    for (const Planet& planet : trappist_system.planetList) {
        y[4*i]     = planet.X[0];
        y[4*i + 1] = planet.X[1];
        y[4*i + 2] = planet.V[0];
        y[4*i + 3] = planet.V[1];
        i++;
    }

    // ========================================================
    // Remove centre-of-mass drift
    // ========================================================
    valarray<double> vCOM = trappist_system.get_barycenter_velocity();
    for (int k = 0; k < nbodies; k++) {
        y[4*k + 2] -= vCOM[0];
        y[4*k + 3] -= vCOM[1];
    }

    // ========================================================
    // Open output file
    // ========================================================
    ofstream fout("trappist.dat");
    if (!fout) {
        cerr << "Error: could not open trappist.dat" << endl;
        return 1;
    }

    fout << fixed << setprecision(10);

    auto writeRow = [&]() {
        fout << t;
        for (int k = 0; k < 4 * nbodies; k++) {
            fout << " " << y[k];
        }
        fout << "\n";
    };

    writeRow();

    // ========================================================
    // Evolve the system
    // ========================================================
    for (int step = 0; step < N; step++) {
        y = RungeKutta4(t, y, h, t + h);
        t += h;
        writeRow();
    }

    fout.close();

    // ========================================================
    // Update system object with final state
    // ========================================================
    trappist_system.star_X[0] = y[0];
    trappist_system.star_X[1] = y[1];
    trappist_system.star_V[0] = y[2];
    trappist_system.star_V[1] = y[3];

    auto it = trappist_system.planetList.begin();
    for (int j = 1; it != trappist_system.planetList.end(); ++it, ++j) {
        it->X[0] = y[4*j];
        it->X[1] = y[4*j + 1];
        it->V[0] = y[4*j + 2];
        it->V[1] = y[4*j + 3];
    }

    cout << "Simulation complete." << endl;
    cout << "Total simulated time: " << t << " s" << endl;
    cout << "Output written to trappist.dat" << endl << endl;

    // final system state
    trappist_system.printInfo();
    cout << endl;

    cout << "Final star state:" << endl;
    cout << "Position: (" << trappist_system.star_X[0] << ", "
         << trappist_system.star_X[1] << ")" << endl;
    cout << "Velocity: (" << trappist_system.star_V[0] << ", "
         << trappist_system.star_V[1] << ")" << endl << endl;

    for (const Planet& planet : trappist_system.planetList) {
        planet.printInfo();
        cout << endl;
    }

    return 0;
}