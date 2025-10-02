#include "core/particle.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

// Declaraciones adelantadas para evitar errores de declaración implícita
void viscosity(std::vector<Particle>& particles);
void computeNavierStokes(std::vector<Particle>& particles);
void boundaryInteraction(std::vector<Particle>& particles, double dx);

// Función unificada para calcular todas las aceleraciones
void acceleration(std::vector<Particle>& particles, double g, double dx) {
    viscosity(particles); // Viscosidad artificial
    computeNavierStokes(particles);
    boundaryInteraction(particles, dx);
    // Nota: applyGravity ya está incluida en computeNavierStokes
}

// Constantes para viscosidad
constexpr double ALPHA_PI = 1.0;
constexpr double BETA_PI = 1.0;
constexpr double EPSILON = 1e-6;

// Viscosidad artificial estilo Monaghan-Gingold
void viscosity(std::vector<Particle>& particles) {
    for (size_t i = 0; i < particles.size(); ++i) {
        Particle& pi = particles[i];
        if (pi.type != 0) continue; // Solo fluido
        for (size_t k = 0; k < pi.neighbors.size(); ++k) {
            int j = pi.neighbors[k];
            Particle& pj = particles[j];

            double xij = pi.pos[0] - pj.pos[0];
            double yij = pi.pos[1] - pj.pos[1];
            double vxij = pi.vel[0] - pj.vel[0];
            double vyij = pi.vel[1] - pj.vel[1];

            double vijrij = vxij * xij + vyij * yij;

            if (vijrij < 0.0) {

                double hij = 0.5 * (pi.h + pj.h);
                double phiij = (hij * vijrij) / (xij * xij + yij * yij + EPSILON);
                
                double cij = 0.5 * (pi.soundVel + pj.soundVel);
                double rhoij = 0.5 * (pi.rho + pj.rho);
                double Piij = (-ALPHA_PI * cij * phiij + BETA_PI * phiij * phiij) / (rhoij + EPSILON);
                
                pi.accel[0] -= pj.mass * Piij * pi.dWx[k];
                pi.accel[1] -= pj.mass * Piij * pi.dWy[k];
                
                double vdw = (pi.vel[0] - pj.vel[0]) * pi.dWx[k] + (pi.vel[1] - pj.vel[1]) * pi.dWy[k];
                pi.dinternalE += 0.5 * pj.mass * Piij * vdw;
            }
        }
    }
}

// Interacción con frontera (simplificada)
void boundaryInteraction(std::vector<Particle>& particles, double dx) {
    int n1 = 12, n2 = 4;
    double r0 = dx / 2.0, D = 0.01;
    for (auto& pi : particles) {
        if (pi.type != 0) continue;
        for (size_t k = 0; k < pi.neighbors.size(); ++k) {
            int j = pi.neighbors[k];
            Particle& pj = particles[j];
            if (pj.type == -1) { // frontera
                double xij = pi.pos[0] - pj.pos[0];
                double yij = pi.pos[1] - pj.pos[1];
                double rij = std::sqrt(xij * xij + yij * yij);
                if (rij < r0) {
                    double factor = D * (std::pow(r0 / rij, n1) - std::pow(r0 / rij, n2)) / (rij * rij + EPSILON);
                    pi.accel[0] += factor * xij;
                    pi.accel[1] += factor * yij;
                }
            }
        }
    }
}

// Aceleración gravitacional
void applyGravity(std::vector<Particle>& particles, double g) {
    for (auto& p : particles) if (p.type == 0) p.accel[1] += g;
}

// Compute Navier-Stokes (llama viscosidad y gravedad)
void computeNavierStokes(std::vector<Particle>& particles) {
    // Aceleración por gradiente de presión
    for (auto& pi : particles) {
        if (pi.type != 0) continue; // solo fluido
        pi.accel[0] = 0.0;
        pi.accel[1] = 0.0;
        pi.dinternalE = 0.0;

        for (size_t k = 0; k < pi.neighbors.size(); k++) {
            int j = pi.neighbors[k];
            Particle& pj = particles[j];

            double pij = (pi.pressure / (pi.rho * pi.rho)) +
                         (pj.pressure / (pj.rho * pj.rho));

            pi.accel[0] -= pj.mass * pij * pi.dWx[k];
            pi.accel[1] -= pj.mass * pij * pi.dWy[k];

            double vdw = (pi.vel[0] - pj.vel[0]) * pi.dWx[k] +
                         (pi.vel[1] - pj.vel[1]) * pi.dWy[k];

            pi.dinternalE += 0.5 * pj.mass * pij * vdw;
        }
    }
}
