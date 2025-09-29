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
    const double h_min_limit = 0.3;
    const double h_max_limit = 3.0;

    for (auto& pi : particles) {
        // Solo partículas de fluido
        if (pi.type != 0) continue;

        if (pi.neighbors.empty()) continue;

        // Factor basado en densidad
        double factor = 1.0 - 0.5 * (pi.rho / rho0 - 1.0);
        if (factor < 0.1) factor = 0.1;

        // Nuevo h candidato
        double h_new = pi.h * factor;

        // Limitar h al rango [0.3, 3.0]
        if (h_new < h_min_limit) h_new = h_min_limit;
        if (h_new > h_max_limit) h_new = h_max_limit;

        pi.h = h_new;
    }
}


