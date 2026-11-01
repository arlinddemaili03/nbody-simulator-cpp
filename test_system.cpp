#include <iostream>
#include <valarray>
#include <list>
#include "Planet.hpp"
#include "System.hpp"

int main() {
    std::valarray<double> pos1 = {1.0e7, 0.0};
    std::valarray<double> vel1 = {0.0, 1000.0};
    Planet p1("Planet1", 1.0e24, 1.0e6, pos1, vel1);

    std::valarray<double> pos2 = {-1.0e7, 0.0};
    std::valarray<double> vel2 = {0.0, -1000.0};
    Planet p2("Planet2", 2.0e24, 1.0e6, pos2, vel2);

    std::list<Planet> planets = {p1, p2};

    std::valarray<double> star_pos = {0.0, 0.0};
    std::valarray<double> star_vel = {0.0, 0.0};

    System sys("TestSystem", "Sun",
               2.0e30, 7.0e8,
               5800, 1.0,
               star_pos, star_vel,
               planets);

    sys.printInfo();

    std::valarray<double> bary_pos = sys.get_barycenter_position();
    std::cout << "Barycenter position: ("
              << bary_pos[0] << ", " << bary_pos[1] << ")" << std::endl;

    std::valarray<double> bary_vel = sys.get_barycenter_velocity();
    std::cout << "Barycenter velocity: ("
              << bary_vel[0] << ", " << bary_vel[1] << ")" << std::endl;

    std::valarray<double> point = {5.0e7, 0.0};
    std::valarray<double> g = sys.gravity_at_point(point);
    std::cout << "Gravity at point: ("
              << g[0] << ", " << g[1] << ")" << std::endl;

    return 0;
}