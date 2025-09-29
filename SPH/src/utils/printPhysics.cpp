#include "utils/printPhysics.h"
#include <iostream>
#include <iomanip>
#include <algorithm>

void printDensityPressure(const std::vector<Particle>& particles, int nPrint) {
    std::cout << "\n=== Debug: Densidad y Presión ===\n";
    for (int i = 0; i < std::min(nPrint, (int)particles.size()); i++) {
        const auto& p = particles[i];
        std::cout << "ID " << p.id
                  << " rho=" << std::fixed << std::setprecision(4) << p.rho
                  << " p="   << std::fixed << std::setprecision(4) << p.pressure
                  << " h="   << std::fixed << std::setprecision(4) << p.h
                  << "\n";
    }
}

void printAccelerations(const std::vector<Particle>& particles, int nPrint) {
    std::cout << "\n=== Debug: Aceleraciones ===\n";
    for (int i = 0; i < std::min(nPrint, (int)particles.size()); i++) {
        const auto& p = particles[i];
        std::cout << "ID " << p.id
                  << " ax=" << std::fixed << std::setprecision(6) << p.accel[0]
                  << " ay=" << std::fixed << std::setprecision(6) << p.accel[1]
                  << " du=" << std::fixed << std::setprecision(6) << p.dinternalE
                  << "\n";
    }
}

#include "utils/printPhysics.h"
#include <iostream>
#include <algorithm> // std::min

void verifyFirstFluidParticle(const std::vector<Particle>& particles) {
    int fluidIndex = -1;

    // Buscar primera partícula fluida
    for (size_t i = 0; i < particles.size(); ++i) {
        if (particles[i].type == 0) { // 0 = fluido
            fluidIndex = static_cast<int>(i);
            break;
        }
    }

    if (fluidIndex >= 0) {
        const Particle& p = particles[fluidIndex];
        std::cout << "=== Verificación partícula de fluido (index=" << fluidIndex << ") ===\n";
        std::cout << "Pos[0]=" << p.pos[0] << ", Pos[1]=" << p.pos[1]
                  << ", H=" << p.h
                  << ", Vel[0]=" << p.vel[0] << ", Vel[1]=" << p.vel[1] << "\n";

        if (!p.neighbors.empty()) {
            std::cout << "Número de vecinos: " << p.neighbors.size() << "\n";
            for (size_t k = 0; k < std::min(p.neighbors.size(), size_t(5)); ++k) {
                int j = p.neighbors[k];
                std::cout << "Vecino " << j
                          << ": dx=" << p.dx[k] << ", dy=" << p.dy[k]
                          << ", r=" << p.r[k]
                          << ", W=" << p.W[k]
                          << ", dWx=" << p.dWx[k]
                          << ", dWy=" << p.dWy[k] << "\n";
            }
        } else {
            std::cout << "No hay vecinos calculados todavía.\n";
        }
    } else {
        std::cout << "No hay partículas de fluido disponibles para verificar.\n";
    }
}

