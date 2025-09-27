#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <limits>
#include "../include/utils.h"

void computeBoundingBox(const std::vector<Particle>& particles,
                        double& xmin, double& xmax,
                        double& ymin, double& ymax) {
    xmin = std::numeric_limits<double>::max();
    xmax = std::numeric_limits<double>::lowest();
    ymin = std::numeric_limits<double>::max();
    ymax = std::numeric_limits<double>::lowest();

    for (const auto& p : particles) {
        if (p.pos[0] < xmin) xmin = p.pos[0];
        if (p.pos[0] > xmax) xmax = p.pos[0];
        if (p.pos[1] < ymin) ymin = p.pos[1];
        if (p.pos[1] > ymax) ymax = p.pos[1];
    }
}

void test_NN(const std::vector<Particle>& particles, int nTests) {
    std::ofstream fTestNN("NN_test.output");
    if (!fTestNN) {
        std::cerr << "No se pudo abrir NN_test.output para escritura\n";
        return;
    }

    srand(time(NULL));
    int nFluid = 0;
    for (const auto& p : particles) {
        if (p.type == 0) nFluid++;
    }

    if (nFluid == 0) {
        std::cerr << "No hay partículas fluidas para test_NN\n";
        return;
    }

    std::cout << "\n=== Test de Vecinos (NN) ===\n";
    for (int k = 0; k < nTests; k++) {
        // Elegir una partícula fluida aleatoria
        int idx;
        do {
            idx = rand() % particles.size();
        } while (particles[idx].type != 0);

        const auto& pi = particles[idx];
        std::cout << "Probando partícula " << pi.id
                  << " con " << pi.neighbors.size() << " vecinos\n";

        // Formato con alineación fija (como en tu código original en C)
        fTestNN << std::setw(16) << pi.id
                << std::setw(16) << std::fixed << std::setprecision(10) << pi.pos[0]
                << std::setw(16) << std::fixed << std::setprecision(10) << pi.pos[1]
                << "\n";

        for (size_t j = 0; j < pi.neighbors.size(); j++) {
            int nj = pi.neighbors[j];
            const auto& pj = particles[nj];
            fTestNN << std::setw(16) << pj.id
                    << std::setw(16) << std::fixed << std::setprecision(10) << pj.pos[0]
                    << std::setw(16) << std::fixed << std::setprecision(10) << pj.pos[1]
                    << "\n";
        }
        fTestNN << "\n";
    }

    fTestNN.close();
    std::cout << "Resultados guardados en NN_test.output\n";
}
