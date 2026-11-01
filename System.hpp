#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include <iostream>
#include <string>
#include <valarray>
#include <list>
#include <cmath>
#include "Planet.hpp"

class System {
public:
    std::string name;
    std::string star_name;
    double star_mass;
    double star_radius;
    double star_temp;
    double star_luminosity;
    std::valarray<double> star_X;
    std::valarray<double> star_V;
    std::list<Planet> planetList;

    // Constructor to initialize all system properties
    System(std::string sname, std::string sstar_name,
           double smass, double sradius,
           double stemp, double slum,
           std::valarray<double> sX,
           std::valarray<double> sV,
           std::list<Planet> plist = {})
        : name(sname), star_name(sstar_name),
          star_mass(smass), star_radius(sradius),
          star_temp(stemp), star_luminosity(slum),
          star_X(sX), star_V(sV), planetList(plist) {}

    // Method to print system information
    void printInfo() {
        std::cout << "System: " << name << std::endl;
        std::cout << "Star: " << star_name << std::endl;
        std::cout << "Mass: " << star_mass << std::endl;
        std::cout << "Planets: " << planetList.size() << std::endl;
    }

    // Method to calculate the barycenter position of the system
    std::valarray<double> get_barycenter_position() {
        std::valarray<double> sum = star_mass * star_X;
        double total_mass = star_mass;

        for (const auto& p : planetList) {
            sum += p.mass * p.X;
            total_mass += p.mass;
        }

        return sum / total_mass;
    }

    // Method to calculate the barycenter velocity of the system
    std::valarray<double> get_barycenter_velocity() {
        std::valarray<double> sum = star_mass * star_V;
        double total_mass = star_mass;

        for (const auto& p : planetList) {
            sum += p.mass * p.V;
            total_mass += p.mass;
        }

        return sum / total_mass;
    }

    // Compute the total gravitational field at point x from the star and all planets
    std::valarray<double> gravity_at_point(std::valarray<double> x) {
        const double G = 6.6743e-11;
        std::valarray<double> g = {0.0, 0.0};

        // Adds the gravitational field due the star
        std::valarray<double> d = star_X - x;
        double r = std::sqrt((d * d).sum());
        g += (G * star_mass / std::pow(r, 3)) * d;

        // Adds the gravitational field due to each planet
        for (const auto& p : planetList) {
            g += p.gravity_at_point(x);
        }

        return g;
    }
};

#endif