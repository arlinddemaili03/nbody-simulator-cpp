#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cassert>
#include <iostream>
#include <valarray>
#include <vector>
#include <cmath>

using namespace std;

/* Physical constants and TRAPPIST-1 system data
   (values supplied with the Queen Mary MTH6150 module materials) */
 
// Earth's radius in meters
#define R_EARTH (6371000.0L)
// Earth's mass in kg
#define M_EARTH (5.9722E+24)
// Solar mass in kg
#define M_SUN (1.989E+30)
// Solar radius in m
#define R_SUN (6.957E+8)
// Solar temperature in K
#define T_SUN (5775.0)
// Newton's constant in SI units (m^3/kg/s^2)
#define G_SI (6.67430E-11)
// Astronomical Units (AU) in meters
#define AU_IN_M (149597870700.0)
// Jupiter's mass in kg
#define M_JUPITER (1.89813E27)
// Jupiter's orbital radius in m
#define D_JUPITER (7.785E11)

const double Y0_EARTH = AU_IN_M;
const double VX0_EARTH = 29789.1;
const double X0_JUP = 0.0;
const double Y0_JUP = -D_JUPITER;
const double VX0_JUP = -13058.4;
const double VY0_JUP = 0.0;

// Gas giant density profile data (used by the interpolation/integration examples)
const double R_GAS_GIANT = 7.0E7; // The gas giant's radius in m
const double r_arr_frac[] = {0.0, 0.125, 0.250, 0.375, 0.5, 0.625, 0.750, 0.875, 1.0};
const double rho_arr_SI[] = {13200, 10200, 7000, 4800, 3000, 1800, 1000, 420, 130};

// Habitable-zone (Kopparapu) coefficients
const double d_inner_coeffs[] = {1.107, 1.332E-4, 1.58E-8, 8.308E-12, 1.931E-15};
const double d_outer_coeffs[] = {0.356, 6.171E-5, 1.698E-9, -3.0045E-12, -5.2983E-16};

// TRAPPIST-1 data
/* Planet orbital radii [AU] */
const double OR_1b_AU = 0.01150;
const double OR_1c_AU = 0.01576;
const double OR_1d_AU = 0.02219;
const double OR_1e_AU = 0.02916;
const double OR_1f_AU = 0.03836;
const double OR_1g_AU = 0.0467;
const double OR_1h_AU = 0.0671;

const double dOR_1b_AU = 0.00028;
const double dOR_1c_AU = 0.00038;
const double dOR_1d_AU = 0.00053;
const double dOR_1e_AU = 0.00070;
const double dOR_1f_AU = 0.00105;
const double dOR_1g_AU = 0.0017;
const double dOR_1h_AU = 0.0019;

/* Planet Masses [kg] */
const double M_1b = 1.02 * M_EARTH;
const double M_1c = 1.31 * M_EARTH;
const double M_1d = 0.39 * M_EARTH;
const double M_1e = 0.69 * M_EARTH;
const double M_1f = 1.04 * M_EARTH;
const double M_1g = 1.32 * M_EARTH;
const double M_1h = 0.33 * M_EARTH;

const double dM_1b = 0.15 * M_EARTH;
const double dM_1c = 0.06 * M_EARTH;
const double dM_1d = 0.02 * M_EARTH;
const double dM_1e = 0.05 * M_EARTH;
const double dM_1f = 0.06 * M_EARTH;
const double dM_1g = 0.04 * M_EARTH;
const double dM_1h = 0.15 * M_EARTH;

/* Star properties of TRAPPIST-1a */
const double M_1a = 0.089 * M_SUN;    // mass in kg
const double dM_1a = 0.006 * M_SUN;   // mass error in kg
const double R_1a = 0.121 * R_SUN;    // stellar radius in m
const double dR_1a = 0.003 * R_SUN;   // stallar radius error in m
const double T_1a = 2511;             // temperature in K
const double dT_1a = 37;              // temperature in K
const double L_1a = 5.22E-4;          // luminosity in units of solar luminosity 
const double dL_1a = 0.19E-4;         // luminosity in units of solar luminosity 


/* Function prototypes */

/* Habitable-zone helpers */
double d_HZ_inner(double starLum, double starTemp, double planetMass);
double d_HZ_outer(double starLum, double starTemp, double planetMass);
bool is_habitable(double starLum, double starTemp, double planetMass, double orbitalRadius);

/* ODE solver */
valarray<double> rhs(double t, const valarray<double> yvec);
valarray<double> RungeKutta4(double t0, valarray<double> y0, double h, double tf);

#endif // CONSTANTS_HPP
