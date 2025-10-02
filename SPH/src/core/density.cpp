#include <iostream>
#include "core/density.h"
#include "core/cubicSplineKernel.h"

void computeDensity(std::vector<Particle>& particles) {
    std::cout << "\n=== Dentro de computeDensity ===\n";
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

    std::cout << "\n=== Dentro de updateSmoothingLength ===\n";
    for (auto& pi : particles) {
        if (pi.type != 0) continue;
        if (pi.neighbors.empty()) continue;

        double old_h = pi.h;

        double factor = 1.0 - 0.5 * (pi.rho / rho0 - 1.0);
        if (factor < 0.1) factor = 0.1;

        double h_new = pi.h * factor;
        if (h_new < 0.3) h_new = 0.3;
        if (h_new > 3.0) h_new = 3.0;

        if (std::abs(h_new - old_h) > 1e-12) {
            std::cout << "[updateSmoothingLength] id=" << pi.id
                      << " old_h=" << old_h
                      << " new_h=" << h_new
                      << " rho=" << pi.rho
                      << " factor=" << factor << "\n";
        }

        pi.h = h_new;
    }
}


