#include "physics/boundary.h"
#include <cmath>
#include <iostream>

void boundaryInteraction(std::vector<Particle>& particles, double dx) {
    int n1 = 12, n2 = 4;
    double r0 = dx / 2.0;
    double D  = 0.01;

    for (auto& pi : particles) {
        if (pi.type != 0) continue; // solo fluido

        for (size_t k = 0; k < pi.neighbors.size(); k++) {
            int j = pi.neighbors[k];
            const Particle& pj = particles[j];

            if (pj.type == -1) { // frontera
                double xij = pi.pos[0] - pj.pos[0];
                double yij = pi.pos[1] - pj.pos[1];
                double rij = std::sqrt(xij * xij + yij * yij);

                if (rij < r0) {
                    double coeff = D * (std::pow(r0 / rij, n1) -
                                        std::pow(r0 / rij, n2)) / (rij * rij);

                    pi.accel[0] += coeff * xij;
                    pi.accel[1] += coeff * yij;
                }
            }
        }
    }
    std::cout << "Boundary interaction computed\n";
}
