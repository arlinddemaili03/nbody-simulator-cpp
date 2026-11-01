#ifndef PLANET_HPP
#define PLANET_HPP

#include <iostream>
#include <string>
#include <valarray>
#include <cmath>

class Planet {
public:
    std::string name;
    double mass;
    double radius;
    std::valarray<double> X;
    std::valarray<double> V;

    // Constructor to initialize all planet's properties
    Planet(std::string pname, double pmass, double pradius,
           std::valarray<double> pX, std::valarray<double> pV)
        : name(pname), mass(pmass), radius(pradius), X(pX), V(pV) {}
    
    // Method to print planet information    
    void printInfo() const {
        std::cout << "Planet: " << name << std::endl;
        std::cout << "Mass: " << mass << " kg" << std::endl;
        std::cout << "Radius: " << radius << " m" << std::endl;
        std::cout << "Position: (" << X[0] << ", " << X[1] << ")" << std::endl;
        std::cout << "Velocity: (" << V[0] << ", " << V[1] << ")" << std::endl;
    }

    // Method to calculate distance from the planet's centre to a given point x
    double distance_to_point(std::valarray<double> x) const {
        std::valarray<double> diff = X - x;
        return std::sqrt((diff * diff).sum());
    }

    // Method to calculate the gravitational acceleration at a given point x
    std::valarray<double> gravity_at_point(std::valarray<double> x) const {
        const double G = 6.6743e-11; // Gravitational constant in m^3 kg^-1 s^-2

        std::valarray<double> diff = X - x;
        double r = std::sqrt((diff * diff).sum());

        return (G * mass / std::pow(r, 3)) * diff;
    }
};

#endif