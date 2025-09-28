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
