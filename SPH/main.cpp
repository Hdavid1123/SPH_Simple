#include <iostream>
#include <limits>
#include "io/particle_io.h"
#include "core/linkedList.h"
#include "core/cubicSplineKernel.h"
#include "utils/utils.h"
#include "core/density.h"
#include "utils/printPhysics.h"
#include "physics/integrator.h"
#include "physics/navierStokes.h"
#include "physics/pressure.h"
#include "utils/output.h"

int main(int argc, char* argv[]) {
    std::string filename = (argc > 1) ? argv[1] : "../data/archivo.txt";

    try {
        auto particles = readParticlesFromFile(filename);
        std::cout << "Se leyeron " << particles.size() << " partículas.\n";

        validateParticles(particles);

        double h_ref = particles[0].h;
        double kappa = 2.0;
        double rho0  = 1000.0;
        int nSteps   = 100;
        double c     = 0.01;
        const double g = -9.81;  // gravedad

        // === Fase 1: malla inicial ===
        auto cells = rebuildGridAndNeighbors(particles, h_ref, kappa);

        computeDensity(particles);
        computePressureEOS(particles, c);
        printDensityPressure(particles, 10);
        test_NN(particles, 20, "NN_before_h_update.output");

        // === Fase 2: actualizar h y vecinos ===
        updateSmoothingLength(particles, rho0);
        cells = rebuildGridAndNeighbors(particles, h_ref, kappa);

        computeDensity(particles);
        computePressureEOS(particles, c);
        printDensityPressure(particles, 10);
        test_NN(particles, 20, "NN_after_h_update.output");

        // === Fase 3: loop de integración ===
        std::cout << "Iniciando integración con " << nSteps << " pasos...\n";

        // Encontrar h mínimo para usar en el cálculo de dt
        double h_min = std::numeric_limits<double>::max();
        for (const auto& p : particles) {
            if (p.h < h_min) h_min = p.h;
        }

        const double nu = 1e-6; // viscosidad cinemática agua (m^2/s)
        // double dt = 0.125 * (h_min * h_min) / nu;
        double dt = 5e-3; // Valor fijo para pruebas iniciales

        std::cout << "dt calculado (Liu 2003) = " << dt 
                  << " usando h_min=" << h_min << "\n";

        int checkInterval = std::max(1, nSteps / 20); // 5% de los pasos

        // ========== Inicio Bucle Temporal ===========
        for (int step = 0; step < nSteps; ++step) {
            std::cout << "Paso " << step << "\n";

            // 1️⃣ Drift (solo partículas de fluido)
            drift(particles, dt);

            // 2️⃣ Actualizar longitud de suavizado y reconstruir malla
            updateSmoothingLength(particles, rho0);
            cells = rebuildGridAndNeighbors(particles, h_ref, kappa);

            // 3️⃣ Recalcular densidad y presión (EOS lineal)
            computeDensity(particles);
            computePressureEOS(particles, c);

            // 4️⃣ Calcular aceleraciones internas (Navier-Stokes)
            computeNavierStokes(particles);

            // 5️⃣ Sumar gravedad a las partículas de fluido
            for (auto& p : particles) {
                if (p.type == 0)
                    p.accel[1] += g;
            }

            // 6️⃣ Kick (solo partículas de fluido)
            kick(particles, dt);

            // 7️⃣ Verificación periódica de una partícula fluida
            if (step % checkInterval == 0) {
                verifyFirstFluidParticle(particles);
            }

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
