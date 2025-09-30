#include <iostream>
#include <limits>
#include <iomanip>  // para std::setprecision
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
#include "physics/boundary.h"

int main(int argc, char* argv[]) {
    std::string filename = (argc > 1) ? argv[1] : "../data/archivo.txt";

    try {
        auto particles = readParticlesFromFile(filename);
        std::cout << "Se leyeron " << particles.size() << " partículas.\n";
        validateParticles(particles);

        double kappa = 2.0;
        double rho0  = 1000.0;
        int nSteps   = 1000;
        double c     = 0.01;
        const double g = -0.987;  // gravedad
        double dt = 1e-3;        // Valor fijo para pruebas iniciales

        // 🔹 Lambda para imprimir cambios en h
        auto log_h_changes = [](const std::vector<Particle>& before,
                                const std::vector<Particle>& after,
                                const std::string& tag) {
            for (size_t i = 0; i < before.size() && i < after.size(); ++i) {
                if (std::abs(before[i].h - after[i].h) > 1e-12) {
                    std::cout << "[h-changed][" << tag << "] id=" << after[i].id
                              << " h_before=" << std::fixed << std::setprecision(6) << before[i].h
                              << " -> h_after=" << std::fixed << std::setprecision(6) << after[i].h
                              << "\n";
                }
            }
        };

        // === Fase 1: inicializar h_min y malla ===
        double h_deb = 0.01; // H fijo para debug

        // Prueba del kernel cúbico
        testKernel(1.0, TWO_D);

        auto before_init = particles;
        auto cells = initializePhase(
            particles, h_deb, kappa, c,
            "Inicial", "NN_before_h_update.output"
        );
        log_h_changes(before_init, particles, "after_initializePhase");

        // === Fase 2: actualizar h y vecinos ===
        auto before_h = particles;
        cells = initializePhase(
            particles, h_deb, kappa, c,
            "Post-h Update", "NN_after_h_update.output"
        );
        log_h_changes(before_h, particles, "after_h_update");

        // === Fase 3: loop de integración ===
        std::cout << "Iniciando integración con " << nSteps << " pasos...\n";

        int checkInterval = std::max(1, nSteps / 20); // 5% de los pasos

        for (int step = 0; step < nSteps; ++step) {
            std::cout << "Paso " << step << "\n";

            // 1️⃣ Drift
            auto before_drift = particles;
            drift(particles, dt);
            log_h_changes(before_drift, particles, "after_drift");

            // 2️⃣ Reconstrucción de malla
            auto before_rebuild = particles;
            cells = rebuildGridAndNeighbors(particles, h_deb, kappa);
            log_h_changes(before_rebuild, particles, "after_rebuildGridAndNeighbors");

            // 3️⃣ Densidad y presión
            auto before_density = particles;
            computeDensity(particles);
            computePressureEOS(particles, c);
            log_h_changes(before_density, particles, "after_density_pressure");

            // 4️⃣ Aceleraciones internas
            auto before_ns = particles;
            computeNavierStokes(particles);
            log_h_changes(before_ns, particles, "after_navier_stokes");

            // 5️⃣ Interacción con frontera
            auto before_boundary = particles;
            boundaryInteraction(particles);
            log_h_changes(before_boundary, particles, "after_boundary");

            // 6️⃣ Gravedad
            auto before_gravity = particles;
            for (auto& p : particles) if (p.type == 0) p.accel[1] += g;
            log_h_changes(before_gravity, particles, "after_gravity");

            // 7️⃣ Kick
            auto before_kick = particles;
            kick(particles, dt);
            log_h_changes(before_kick, particles, "after_kick");

            // 8️⃣ Verificación periódica
            if (step % checkInterval == 0) {
                auto before_verify = particles;
                verifyFirstFluidParticle(particles);
                log_h_changes(before_verify, particles, "after_verify");
            }

            // 9️⃣ Guardar estado
            auto before_print = particles;
            printState(particles, step);
            log_h_changes(before_print, particles, "after_printState");
        }

        std::cout << "Integración finalizada.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
