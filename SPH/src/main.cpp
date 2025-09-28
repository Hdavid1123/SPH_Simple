#include <iostream>
#include "../include/particle_io.h"
#include "../include/linkedList.h"
#include "../include/cubicSplineKernel.h"
#include "../include/utils.h"
#include "../include/density.h"   // nueva cabecera para densidad

int main(int argc, char* argv[]) {
    std::string filename = (argc > 1) ? argv[1] : "../data/archivo.txt";

    try {
        auto particles = readParticlesFromFile(filename);
        std::cout << "Se leyeron " << particles.size() << " partículas.\n";

        validateParticles(particles);

        double xmin, xmax, ymin, ymax;
        computeBoundingBox(particles, xmin, xmax, ymin, ymax);

        double h_ref = particles[0].h;
        double kappa = 2.0;
        double rho0  = 1000.0;  // densidad de referencia (agua)

        // === Fase 1: con h inicial ===
        auto cells = buildCellGrid(xmin - h_ref, xmax + h_ref,
                                   ymin - h_ref, ymax + h_ref,
                                   h_ref, kappa);

        assignParticlesToCells(cells, particles, h_ref, kappa);
        findNeighbors(cells, particles, kappa);

        computeDensity(particles);
        test_NN(particles, 20, "NN_before_h_update.output");

        // === Fase 2: actualizar h y recalcular vecinos ===
        updateSmoothingLength(particles, rho0);

        cells = buildCellGrid(xmin - h_ref, xmax + h_ref,
                              ymin - h_ref, ymax + h_ref,
                              h_ref, kappa);

        assignParticlesToCells(cells, particles, h_ref, kappa);
        findNeighbors(cells, particles, kappa);

        computeDensity(particles);
        test_NN(particles, 20, "NN_after_h_update.output");

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
