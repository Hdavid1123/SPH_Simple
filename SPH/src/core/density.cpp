#include "core/density.h"
#include "core/cubicSplineKernel.h"

void computeDensity(std::vector<Particle>& particles) {
    for (auto& pi : particles) {
        double rho = 0.0;

        // self contribution
        rho += pi.mass * cubicSplineKernel(0.0, pi.h, TWO_D);

        // vecinos
        for (size_t n = 0; n < pi.neighbors.size(); n++) {
            int j = pi.neighbors[n];
            rho += particles[j].mass * pi.W[n];
        }

        pi.rho = rho;
    }
}

void updateSmoothingLength(std::vector<Particle>& particles, double rho0) {
    for (auto& pi : particles) {
        double factor = 1.0 - 0.5 * (pi.rho / rho0 - 1.0);

        // proteger contra valores negativos o muy pequeños
        if (factor < 0.1) factor = 0.1;
        pi.h *= factor;
    }
}
