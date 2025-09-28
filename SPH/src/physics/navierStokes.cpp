#include "physics/navierStokes.h"
#include "physics/pressure.h"   // usa computePressureEOS o computePressureTait
#include <vector>

// Versión simple, equivalente a tu código viejo pero con std::vector
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
