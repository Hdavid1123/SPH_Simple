#include "physics/boundary.h"
#include <cmath>
#include <iostream>

void boundaryInteraction(std::vector<Particle>& particles) {
    int n1 = 12, n2 = 4;
    double D  = 0.01;

    for (auto& pi : particles) {
        if (pi.type != 0) continue; // solo fluido

        double r0_sum = 0.0;
        int n_boundary = 0;

        // Sumar distancias a vecinos frontera (type == 1)
        for (size_t k = 0; k < pi.neighbors.size(); k++) {
            int j = pi.neighbors[k];
            const Particle& pj = particles[j];
            if (pj.type == 1) { // frontera
                r0_sum += pi.r[k];
                n_boundary++;
            }
        }

        if (n_boundary == 0) continue; // sin vecinos de frontera

        //double r0_avg = r0_sum / n_boundary;
        double r0_avg = 0.7 * (pi.h); // usar h en lugar de distancia promedio

        // Aplicar la fuerza repulsiva a los vecinos frontera
        for (size_t k = 0; k < pi.neighbors.size(); k++) {
            int j = pi.neighbors[k];
            const Particle& pj = particles[j];
            if (pj.type != 1) continue;

            double rij = pi.r[k];
            if (rij < r0_avg) {
                double coeff = D * (std::pow(r0_avg / rij, n1) - std::pow(r0_avg / rij, n2)) / (rij * rij);
                pi.accel[0] += coeff * pi.dx[k];
                pi.accel[1] += coeff * pi.dy[k];
            }
        }
    }

    std::cout << "Boundary interaction computed\n";
}

void boundaryInteractionCuadratic(std::vector<Particle>& particles) {
    double ka = 1e4;
    double d = 0.008;
    double n = 2.0;

    for (auto& pi : particles) {
        if (pi.type != 0) continue; // solo fluido

        double r0_sum = 0.0;
        int n_boundary = 0;

        // Sumar distancias a vecinos frontera (type == 1)
        for (size_t k = 0; k < pi.neighbors.size(); k++) {
            int j = pi.neighbors[k];
            const Particle& pj = particles[j];
            if (pj.type == 1) { // frontera
                r0_sum += pi.r[k];
                n_boundary++;
            }
        }

        if (n_boundary == 0) continue; // sin vecinos de frontera

        // Aplicar la fuerza repulsiva a los vecinos frontera
        for (size_t k = 0; k < pi.neighbors.size(); k++) {
            int j = pi.neighbors[k];
            const Particle& pj = particles[j];
            if (pj.type != 1) continue;

            double rij = pi.r[k];
            std::cout<<rij<<d<<"\n";
            if (rij < d && rij > 1e-12) { // evitar división por cero
                double coeff = - ka * (std::pow(d - rij, n)) / (rij);
                pi.accel[0] += coeff * pi.dx[k];
                pi.accel[1] += coeff * pi.dy[k];
            }
        }
    }

    std::cout << "Boundary interaction computed\n";
}


