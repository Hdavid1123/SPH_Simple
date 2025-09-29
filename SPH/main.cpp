#include <iostream>
#include <limits>
#include "io/particle_io.h"
#include "core/linkedList.h"
#include "core/cubicSplineKernel.h"
#include "core/density.h"
#include "physics/integrator.h"
#include "physics/navierStokes.h"
#include "physics/pressure.h"
#include "utils/output.h"
#include "utils/printPhysics.h"
#include "utils/simulationUtils.h"

int main(int argc, char* argv[]) {
    std::string filename = (argc > 1) ? argv[1] : "../data/archivo.txt";

    try {
        auto particles = readParticlesFromFile(filename);
        std::cout << "Se leyeron " << particles.size() << " partículas.\n";
        validateParticles(particles);

        double kappa = 2.0;
        double rho0  = 1000.0;
        int nSteps   = 100;
        double c     = 0.01;
        const double g = -9.81;  // gravedad
        double dt = 5e-3;        // Valor fijo para pruebas iniciales

        // === Fase 1: inicializar h_min y malla ===
        double h_min = std::numeric_limits<double>::max();
        for (const auto& p : particles) if (p.type == 0) h_min = std::min(h_min, p.h);

        auto cells = initializePhase(
            particles, h_min, kappa, c,
            "Inicial", "NN_before_h_update.output"
        );

        // === Fase 2: actualizar h y vecinos ===
        updateSmoothingLength(particles, rho0);

        // recalcular h_min después de actualizar h
        h_min = std::numeric_limits<double>::max();
        for (const auto& p : particles) if (p.type == 0) h_min = std::min(h_min, p.h);

        cells = initializePhase(
            particles, h_min, kappa, c,
            "Post-h Update", "NN_after_h_update.output"
        );

        // === Fase 3: loop de integración ===
        std::cout << "Iniciando integración con " << nSteps << " pasos...\n";

        int checkInterval = std::max(1, nSteps / 20); // 5% de los pasos

        for (int step = 0; step < nSteps; ++step) {
            std::cout << "Paso " << step << "\n";

            // 1️⃣ Drift (solo partículas de fluido)
            drift(particles, dt);

            // 2️⃣ Actualizar longitud de suavizado y reconstruir malla
            updateSmoothingLength(particles, rho0);

            // recalcular h_min para la malla
            h_min = std::numeric_limits<double>::max();
            for (const auto& p : particles) if (p.type == 0) h_min = std::min(h_min, p.h);

            cells = rebuildGridAndNeighbors(particles, h_min, kappa);

            // 3️⃣ Recalcular densidad y presión
            computeDensity(particles);
            computePressureEOS(particles, c);

            // 4️⃣ Aceleraciones internas (Navier-Stokes)
            computeNavierStokes(particles);

            // 5️⃣ Sumar gravedad a partículas de fluido
            for (auto& p : particles) if (p.type == 0) p.accel[1] += g;

            // 6️⃣ Kick
            kick(particles, dt);

            // 7️⃣ Verificación periódica
            if (step % checkInterval == 0) verifyFirstFluidParticle(particles);

            // 8️⃣ Guardar estado
            printState(particles, step);
        }

        std::cout << "Integración finalizada.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
