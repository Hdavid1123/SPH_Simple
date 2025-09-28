#include "physics/pressure.h"
#include <cmath>

// EOS lineal: p = c^2 * rho
void computePressureEOS(std::vector<Particle>& particles, double c) {
    for (auto& pi : particles) {
        pi.soundVel = c;  // velocidad del sonido local (constante aquí)
        pi.pressure = c * c * pi.rho;
    }
}

// Tait EOS: p = c0^2 * rho0 * [ (rho / rho0)^gamma - 1 ]
void computePressureTait(std::vector<Particle>& particles,
                         double rho0, double c0, double gamma) {
    for (auto& pi : particles) {
        pi.soundVel = c0;  // velocidad del sonido artificial global
        pi.pressure = c0 * c0 * rho0 * (std::pow(pi.rho / rho0, gamma) - 1.0);
    }
}
