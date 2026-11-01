#include <iostream>
#include <iomanip>
#include <fstream>
#include <valarray>
#include <cmath>
#include <cstdlib>
#include "constants.hpp"

using namespace std;

// ============================================================
// 2-body RHS (kept here so we can still call it if needed)
// State vector: [xs, ys, vxs, vys, xe, ye, vxe, vye]
// ============================================================
valarray<double> rhs2Body(double t, const valarray<double> yvec) {
    (void)t;

    valarray<double> dydt(8);

    double xs  = yvec[0];
    double ys  = yvec[1];
    double vxs = yvec[2];
    double vys = yvec[3];

    double xe  = yvec[4];
    double ye  = yvec[5];
    double vxe = yvec[6];
    double vye = yvec[7];

    double dx = xe - xs;
    double dy = ye - ys;
    double r  = sqrt(dx*dx + dy*dy);
    double r3 = r*r*r;

    // Position derivatives
    dydt[0] = vxs;
    dydt[1] = vys;
    dydt[4] = vxe;
    dydt[5] = vye;

    // Velocity derivatives (gravity)
    dydt[2] =  G_SI * M_EARTH * dx / r3;
    dydt[3] =  G_SI * M_EARTH * dy / r3;
    dydt[6] = -G_SI * M_SUN   * dx / r3;
    dydt[7] = -G_SI * M_SUN   * dy / r3;

    return dydt;
}

// ============================================================
// 3-body RHS for Sun-Earth-Jupiter
// State vector: [xs, ys, vxs, vys, xe, ye, vxe, vye, xj, yj, vxj, vyj]
// ============================================================
valarray<double> rhs_3Body(double t, const valarray<double> yvec) {
    (void)t;

    valarray<double> dydt(12);

    double xs  = yvec[0],  ys  = yvec[1];
    double vxs = yvec[2],  vys = yvec[3];
    double xe  = yvec[4],  ye  = yvec[5];
    double vxe = yvec[6],  vye = yvec[7];
    double xj  = yvec[8],  yj  = yvec[9];
    double vxj = yvec[10], vyj = yvec[11];

    // Pairwise separations
    double dx_se = xe - xs, dy_se = ye - ys;
    double r_se  = sqrt(dx_se*dx_se + dy_se*dy_se);
    double r_se3 = r_se*r_se*r_se;

    double dx_sj = xj - xs, dy_sj = yj - ys;
    double r_sj  = sqrt(dx_sj*dx_sj + dy_sj*dy_sj);
    double r_sj3 = r_sj*r_sj*r_sj;

    double dx_ej = xj - xe, dy_ej = yj - ye;
    double r_ej  = sqrt(dx_ej*dx_ej + dy_ej*dy_ej);
    double r_ej3 = r_ej*r_ej*r_ej;

    // Position derivatives
    dydt[0]  = vxs;  dydt[1]  = vys;
    dydt[4]  = vxe;  dydt[5]  = vye;
    dydt[8]  = vxj;  dydt[9]  = vyj;

    // Sun's acceleration: pulled by Earth (+dx_se direction) and Jupiter
    dydt[2] =  G_SI * M_EARTH   * dx_se / r_se3
             + G_SI * M_JUPITER * dx_sj / r_sj3;
    dydt[3] =  G_SI * M_EARTH   * dy_se / r_se3
             + G_SI * M_JUPITER * dy_sj / r_sj3;

    // Earth's acceleration: pulled by Sun (-dx_se) and Jupiter (+dx_ej)
    dydt[6] = -G_SI * M_SUN     * dx_se / r_se3
             + G_SI * M_JUPITER * dx_ej / r_ej3;
    dydt[7] = -G_SI * M_SUN     * dy_se / r_se3
             + G_SI * M_JUPITER * dy_ej / r_ej3;

    // Jupiter's acceleration: pulled by Sun (-dx_sj) and Earth (-dx_ej)
    dydt[10] = -G_SI * M_SUN    * dx_sj / r_sj3
              - G_SI * M_EARTH  * dx_ej / r_ej3;
    dydt[11] = -G_SI * M_SUN    * dy_sj / r_sj3
              - G_SI * M_EARTH  * dy_ej / r_ej3;

    return dydt;
}

// ============================================================
// Wrapper used by RK4 in rk4.cpp.
// Points to the 3-body RHS for the Sun-Earth-Jupiter simulation.
// ============================================================
valarray<double> rhs(double t, const valarray<double> yvec) {
    return rhs_3Body(t, yvec);
}

int main(int argc, char* argv[]) {
    // Command-line args: step size h and number of steps N
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <h> <N>" << endl;
        return 1;
    }

    double h = atof(argv[1]);
    int N    = atoi(argv[2]);
    double t = 0.0;

    // ========================================================
    // 3-body state vector for Sun, Earth, Jupiter
    // Layout: [xs, ys, vxs, vys, xe, ye, vxe, vye, xj, yj, vxj, vyj]
    // ========================================================
    valarray<double> y(12);

    // Sun at origin, at rest
    y[0] = 0.0;   y[1] = 0.0;
    y[2] = 0.0;   y[3] = 0.0;

    // Earth at (0, Y0_EARTH) moving in +x
    y[4] = 0.0;        y[5] = Y0_EARTH;
    y[6] = VX0_EARTH;  y[7] = 0.0;

    // Jupiter at (X0_JUP, Y0_JUP) with given velocity
    y[8]  = X0_JUP;    y[9]  = Y0_JUP;
    y[10] = VX0_JUP;   y[11] = VY0_JUP;

    // ========================================================
    // remove centre-of-mass drift by subtracting
    // the COM velocity from each body's initial velocity.
    // ========================================================
    double total_mass = M_SUN + M_EARTH + M_JUPITER;
    double vcom_x = (M_SUN * y[2] + M_EARTH * y[6] + M_JUPITER * y[10]) / total_mass;
    double vcom_y = (M_SUN * y[3] + M_EARTH * y[7] + M_JUPITER * y[11]) / total_mass;

    y[2]  -= vcom_x;   y[3]  -= vcom_y;
    y[6]  -= vcom_x;   y[7]  -= vcom_y;
    y[10] -= vcom_x;   y[11] -= vcom_y;

    // Save initial state for later comparison
    valarray<double> y_initial = y;

    // ========================================================
    // Output every 10000 seconds (keeps file size manageable)
    // ========================================================
    int output_every = static_cast<int>(10000.0 / h);
    if (output_every < 1) output_every = 1;

    ofstream fout("sun_earth_jupiter.dat");
    if (!fout) {
        cerr << "Error: could not open sun_earth_jupiter.dat" << endl;
        return 1;
    }

    fout << fixed << setprecision(10);

    // Columns: t xs ys vxs vys xe ye vxe vye xj yj vxj vyj
    fout << t << " "
         << y[0]  << " " << y[1]  << " " << y[2]  << " " << y[3]  << " "
         << y[4]  << " " << y[5]  << " " << y[6]  << " " << y[7]  << " "
         << y[8]  << " " << y[9]  << " " << y[10] << " " << y[11] << "\n";

    for (int i = 0; i < N; ++i) {
        y = RungeKutta4(t, y, h, t + h);
        t += h;

        if ((i + 1) % output_every == 0) {
            fout << t << " "
                 << y[0]  << " " << y[1]  << " " << y[2]  << " " << y[3]  << " "
                 << y[4]  << " " << y[5]  << " " << y[6]  << " " << y[7]  << " "
                 << y[8]  << " " << y[9]  << " " << y[10] << " " << y[11] << "\n";
        }
    }

    fout.close();

    // ========================================================
    // Print initial and final states for comparison
    // ========================================================
    cout << scientific << setprecision(6);
    cout << "Simulation complete." << endl;
    cout << "Total simulated time: " << t << " s  (~ "
         << t / 3.15576e7 << " years)" << endl;
    cout << "Output written to sun_earth_jupiter.dat" << endl;

    cout << "\n=== Earth state ===" << endl;
    cout << "Initial: x=" << y_initial[4] << ", y=" << y_initial[5]
         << ", vx=" << y_initial[6] << ", vy=" << y_initial[7] << endl;
    cout << "Final:   x=" << y[4]         << ", y=" << y[5]
         << ", vx=" << y[6]         << ", vy=" << y[7]         << endl;
    cout << "Delta:   dx=" << y[4]-y_initial[4]
         << ", dy=" << y[5]-y_initial[5]
         << ", dvx=" << y[6]-y_initial[6]
         << ", dvy=" << y[7]-y_initial[7] << endl;

    cout << "\n=== Jupiter state ===" << endl;
    cout << "Initial: x=" << y_initial[8]  << ", y=" << y_initial[9]
         << ", vx=" << y_initial[10] << ", vy=" << y_initial[11] << endl;
    cout << "Final:   x=" << y[8]          << ", y=" << y[9]
         << ", vx=" << y[10]         << ", vy=" << y[11]         << endl;

    cout << "\n=== Sun state ===" << endl;
    cout << "Initial: x=" << y_initial[0]  << ", y=" << y_initial[1]
         << ", vx=" << y_initial[2]  << ", vy=" << y_initial[3] << endl;
    cout << "Final:   x=" << y[0]          << ", y=" << y[1]
         << ", vx=" << y[2]          << ", vy=" << y[3] << endl;

    return 0;
}