# N-Body Planetary Simulator (C++)

A small C++17 project that models planetary systems with Newtonian gravity. It integrates the equations of motion with a fourth-order Runge–Kutta (RK4) solver, simulates the seven-planet TRAPPIST-1 system and a Sun–Earth–Jupiter system, and uses a Monte Carlo simulation to estimate how likely each TRAPPIST-1 planet is to sit in its star's habitable zone once measurement uncertainty is taken into account.

I built this during my BSc Mathematics at Queen Mary University of London, as part of the *Numerical Computing with C and C++* module.

![TRAPPIST-1 orbits](plots/Trappist%20Graph%20plotted.png)

## What's in here

| File | What it does |
|---|---|
| `Planet.hpp` | `Planet` class: name, mass, radius, position and velocity vectors, distance to a point and the gravitational field it produces there |
| `System.hpp` | `System` class: a star plus a `std::list<Planet>`, with barycentre position/velocity and the total gravitational field at a point |
| `rk4.cpp` | Classic RK4 integrator for a system of first-order ODEs `y' = f(t, y)` |
| `trappist.cpp` | General N-body right-hand side; simulates the star and all 7 TRAPPIST-1 planets |
| `sun_earth_jupiter.cpp` | Two- and three-body simulations of the Sun, Earth and Jupiter |
| `monte_carlo_habitability.cpp` | 10 million-sample Monte Carlo estimate of the probability each planet is habitable |
| `test_planet.cpp`, `test_system.cpp` | Small checks for the two classes |
| `constants.hpp` | Physical constants and system parameters |

## Building and running

Requires a C++17 compiler (tested with g++).

```bash
make                          # builds everything

./trappist 100 10000          # step size h = 100 s, N = 10,000 steps -> trappist.dat
./sun_earth_jupiter 100 10000 # -> sun_earth_jupiter.dat
./monte_carlo_habitability    # prints a habitability probability for each planet
./test_planet && ./test_system
```

The simulations write one row per output step (time, then x, y, vx, vy for each body), which can be plotted with Python/matplotlib, Excel or anything else.

## How it works

**State vector.** Every body contributes four values `[x, y, vx, vy]`, so an N-body system is one `std::valarray<double>` of length `4N`. The right-hand side returns the velocities and the accelerations, where each body's acceleration is the sum of `G m_j (r_j − r_i) / |r_j − r_i|^3` over all the other bodies. Because the RHS is written for a general N, the same integrator handles 2, 3 or 8 bodies.

**Removing barycentre drift.** If the star starts at rest while the planets are moving, the system has net momentum and the whole thing slowly drifts across the plot. Before integrating, I subtract the centre-of-mass velocity from every body, so the barycentre stays fixed and the orbits close properly.

**Sanity checks.** `sun_earth_jupiter` prints each body's initial and final state, so you can check that orbits close up after a full period and that the barycentre stays put. Halving the step size and comparing results is a quick way to check RK4's fourth-order convergence.

**Monte Carlo habitability.** The habitable-zone edges come from the Kopparapu et al. fit, which depends on the star's luminosity and temperature and on the planet's mass. All of those, plus each planet's orbital radius, are only known to within measurement error. The program treats the four quantities as independent Gaussian random variables, draws 10 million samples per planet, and counts how often the planet lands between the inner and outer edges of the zone. Results:

| Planet | P(habitable) |
|---|---|
| TRAPPIST-1e | ≈ 1.00 |
| TRAPPIST-1f | ≈ 1.00 |
| TRAPPIST-1g | ≈ 0.93 |
| b, c, d, h | ≈ 0 |

## Credits

The RK4 routine in `rk4.cpp`, and the physical constants and system parameters in `constants.hpp`, are based on material supplied with the Queen Mary module. The classes, N-body right-hand sides, simulations and Monte Carlo analysis are my own work.
