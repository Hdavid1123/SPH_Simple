#include <iostream>
#include <limits>
#include <iomanip>  // para std::setprecision
#include "io/particle_io.h"
#include "core/linkedList.h"
#include "core/cubicSplineKernel.h"
#include "core/density.h"
#include "physics/integrator.h"
#include "physics/pressure.h"
#include "utils/output.h"
#include "utils/printPhysics.h"
#include "utils/simulationUtils.h"
#include "physics/acceleration.h"

int main(int argc, char* argv[]) {
    std::string filename = (argc > 1) ? argv[1] : "../data/archivo.txt";

    try {
        auto particles = readParticlesFromFile(filename);
        std::cout << "Se leyeron " << particles.size() << " partículas.\n";
        //validateParticles(particles);

        double kappa = 2.0;
        double rho0  = 1000.0;
        int nSteps   = 1000;
        double c     = 0.01;
        const double g = -9.87;  // gravedad
        double dt = 1e-3;        // Valor fijo para pruebas iniciales
        double H = 0.15;         // Altura de la columna debe calcularse de mejor manera

        // 🔹 Helper para imprimir los primeros N valores de h
        auto print_h_values = [](const std::vector<Particle>& particles,
                                 const std::string& tag, int n=5) {
            std::cout << "\n--- " << tag << " ---\n";
            for (int i = 0; i < std::min(n, (int)particles.size()); i++) {
                std::cout << "id=" << particles[i].id
                          << " h=" << std::fixed << std::setprecision(6) << particles[i].h
                          << "\n";
            }
        };

        // === Fase 1: inicializar h_min y malla ===
        double h_deb = 0.01; // H fijo para debug

        // Prueba del kernel cúbico
        testKernel(1.0, TWO_D);

        print_h_values(particles, "Después de leer archivo");

        auto before_init = particles;
        auto cells = initializePhase(
            particles, h_deb, kappa, c,
            "Inicial", "NN_before_h_update.output"
        );
        print_h_values(particles, "Después de initializePhase [Inicial]");

        // === Fase 2: actualizar h y vecinos ===
        auto before_h = particles;
        cells = initializePhase(
            particles, h_deb, kappa, c,
            "Post-h Update", "NN_after_h_update.output"
        );
        print_h_values(particles, "Después de initializePhase [Post-h Update]");

        // === Fase 3: loop de integración ===
        std::cout << "Iniciando integración con " << nSteps << " pasos...\n";

        int checkInterval = std::max(1, nSteps / 20); // 5% de los pasos

        for (int step = 0; step < nSteps; ++step) {
            std::cout << "Paso " << step << "\n";
            //print_h_values(particles, "Inicio paso " + std::to_string(step));

            // 1️⃣ Drift
            auto before_drift = particles;
            drift(particles, dt);
            //print_h_values(particles, "Después de Drift");

            // 2️⃣ Reconstrucción de malla
            auto before_rebuild = particles;
            cells = rebuildGridAndNeighbors(particles, h_deb, kappa);
            //print_h_values(particles, "Después de Rebuild");
            std::cout << "Número de celdas: " << cells.size() << "\n";
            
            // 3️⃣ Densidad y presión
            auto before_density = particles;
            computeDensity(particles);
            computePressureKorzani(particles, H, g, rho0);
            //computePressureMonaghan(particles, rho0);
            //computePressureEOS(particles, c)
            //print_h_values(particles, "Después de Density+Pressure");

            // 4️⃣ Aceleraciones internas, frontera y gravedad
            auto before_accel = particles;
            acceleration(particles, g, h_deb);
            //print_h_values(particles, "Después de Acceleration");

            // 7️⃣ Kick
            auto before_kick = particles;
            kick(particles, dt);
            //print_h_values(particles, "Después de Kick");

            // 8️⃣ Verificación periódica
            if (step % checkInterval == 0) {
                auto before_verify = particles;
                verifyFirstFluidParticle(particles);
                print_h_values(particles, "Después de Verify");
            }

            // 9️⃣ Guardar estado
            auto before_print = particles;
            printState(particles, step);
            //print_h_values(particles, "Después de printState");

            // 10️⃣ Liberar memoria temporal
            freeParticleNeighbors(particles);
        }

        std::cout << "Integración finalizada.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
