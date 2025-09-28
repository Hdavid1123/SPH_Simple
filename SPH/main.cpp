#include <iostream>
#include "io/particle_io.h"
#include "core/linkedList.h"
#include "core/cubicSplineKernel.h"
#include "utils/utils.h"
#include "core/density.h"
#include "utils/printPhysics.h"
#include "physics/integrator.h"
#include "physics/navierStokes.h"
#include "physics/pressure.h"

int main(int argc, char* argv[]) {
    std::string filename = (argc > 1) ? argv[1] : "../data/archivo.txt";

    try {
        auto particles = readParticlesFromFile(filename);
        std::cout << "Se leyeron " << particles.size() << " partículas.\n";

        validateParticles(particles);

        double h_ref = particles[0].h;
        double kappa = 2.0;
        double rho0  = 1000.0;
        double dt    = 0.001;
        int nSteps   = 20;
        double c     = 0.01;
        const double g = -9.81;  // gravedad

        // === Fase 1: malla inicial ===
        double xmin, xmax, ymin, ymax;
        computeBoundingBox(particles, xmin, xmax, ymin, ymax);

        auto cells = buildCellGrid(xmin - h_ref, xmax + h_ref,
                                   ymin - h_ref, ymax + h_ref,
                                   h_ref, kappa);

        assignParticlesToCells(cells, particles, h_ref, kappa);
        findNeighbors(cells, particles, kappa);

        computeDensity(particles);
        computePressureEOS(particles, c);
        printDensityPressure(particles, 10);
        test_NN(particles, 20, "NN_before_h_update.output");

        // === Fase 2: actualizar h y vecinos ===
        updateSmoothingLength(particles, rho0);

        computeBoundingBox(particles, xmin, xmax, ymin, ymax);
        cells = buildCellGrid(xmin - h_ref, xmax + h_ref,
                              ymin - h_ref, ymax + h_ref,
                              h_ref, kappa);

        assignParticlesToCells(cells, particles, h_ref, kappa);
        findNeighbors(cells, particles, kappa);

        computeDensity(particles);
        computePressureEOS(particles, c);
        printDensityPressure(particles, 10);
        test_NN(particles, 20, "NN_after_h_update.output");

        // === Fase 3: loop de integración ===
        std::cout << "Iniciando integración con " << nSteps << " pasos...\n";

        for (int step = 0; step < nSteps; ++step) {
            std::cout << "Paso " << step << "\n";

            // 1️⃣ Drift (solo partículas de fluido)
            drift(particles, dt);

            // 2️⃣ Actualizar longitud de suavizado y reconstruir malla
            updateSmoothingLength(particles, rho0);

            computeBoundingBox(particles, xmin, xmax, ymin, ymax);
            cells = buildCellGrid(xmin - h_ref, xmax + h_ref,
                                  ymin - h_ref, ymax + h_ref,
                                  h_ref, kappa);
            assignParticlesToCells(cells, particles, h_ref, kappa);
            findNeighbors(cells, particles, kappa);

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

            // 7️⃣ (Opcional) buscar la primera partícula de tipo fluido
            int fluidIndex = -1;
            for (size_t i = 0; i < particles.size(); ++i) {
                if (particles[i].type == 0) {
                    fluidIndex = i;
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

        std::cout << "Integración finalizada.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
