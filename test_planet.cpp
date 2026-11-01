#include <iostream>
#include <valarray>
#include "Planet.hpp"

int main() {
    std::valarray<double> pos = {0.0, 0.0};
    std::valarray<double> vel = {0.0, 1000.0};

    Planet earth("Earth", 5.972e24, 6.371e6, pos, vel);

    earth.printInfo();

    std::valarray<double> point = {1.0e7, 0.0};

    std::cout << "Distance to point: "
              << earth.distance_to_point(point) << " m" << std::endl;

    std::valarray<double> g = earth.gravity_at_point(point);

    std::cout << "Gravity at point: (" << g[0] << ", " << g[1] << ")" << std::endl;

    return 0;
}