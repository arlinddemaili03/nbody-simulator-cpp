CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra
BINS = trappist sun_earth_jupiter monte_carlo_habitability test_planet test_system

all: $(BINS)

trappist: trappist.cpp rk4.cpp Planet.hpp System.hpp constants.hpp
	$(CXX) $(CXXFLAGS) -o $@ trappist.cpp rk4.cpp

sun_earth_jupiter: sun_earth_jupiter.cpp rk4.cpp constants.hpp
	$(CXX) $(CXXFLAGS) -o $@ sun_earth_jupiter.cpp rk4.cpp

monte_carlo_habitability: monte_carlo_habitability.cpp constants.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

test_planet: test_planet.cpp Planet.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

test_system: test_system.cpp System.hpp Planet.hpp
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(BINS) *.dat

.PHONY: all clean
