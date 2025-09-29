#include <iomanip>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <filesystem>
#include "utils/utils.h"

void test_NN(const std::vector<Particle>& particles, int nTests,
             const std::string& filename) {
    namespace fs = std::filesystem;

    // Crear carpeta test_results si no existe
    fs::create_directories("test_results");

    // Ruta completa dentro de test_results
    std::string filepath = "test_results/" + filename;

    std::ofstream fTestNN(filepath);
    if (!fTestNN) {
        std::cerr << "No se pudo abrir " << filepath << " para escritura\n";
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
        int idx;
        do {
            idx = rand() % particles.size();
        } while (particles[idx].type != 0);

        const auto& pi = particles[idx];
        std::cout << "Probando partícula " << pi.id
                  << " con " << pi.neighbors.size() << " vecinos\n";

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
    std::cout << "Resultados guardados en " << filepath << "\n";
}