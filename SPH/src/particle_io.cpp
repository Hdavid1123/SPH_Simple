#include "../include/particle_io.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <limits>

std::vector<Particle> readParticlesFromFile(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Error: no se pudo abrir " + filename);
    }

    std::vector<Particle> particles;
    std::string line;

    // Saltar la primera línea (cabecera)
    std::getline(infile, line);

    while (std::getline(infile, line)) {
        if (line.empty()) continue; // saltar líneas vacías
        std::istringstream iss(line);

        Particle p;
        double posx, posy, dx, dy;

        if (!(iss >> p.id >> posx >> posy >> p.h >> p.type >> p.mass >> dx >> dy)) {
            std::cerr << "Error en formato de línea: " << line << "\n";
            continue;
        }

        p.pos = {posx, posy};
        p.vel = {0.0, 0.0};
        p.accel = {0.0, 0.0};
        p.rho = 0.0;
        p.pressure = 0.0;
        p.soundVel = 0.0;
        p.internalE = 0.0;
        p.dinternalE = 0.0;

        particles.push_back(p);
    }

    return particles;
}

void validateParticles(const std::vector<Particle>& particles) {
    if (particles.empty()) {
        std::cerr << "Advertencia: no se leyeron partículas.\n";
        return;
    }

    // Rango de valores para reporte
    double min_h = std::numeric_limits<double>::max();
    double max_h = 0.0;
    double min_mass = std::numeric_limits<double>::max();
    double max_mass = 0.0;

    int invalid_count = 0;

    for (const auto& p : particles) {
        if (p.h <= 0 || p.mass <= 0) {
            invalid_count++;
        }

        if (p.h < min_h) min_h = p.h;
        if (p.h > max_h) max_h = p.h;
        if (p.mass < min_mass) min_mass = p.mass;
        if (p.mass > max_mass) max_mass = p.mass;
    }

    std::cout << "Resumen de validación:\n";
    std::cout << " - Partículas leídas: " << particles.size() << "\n";
    std::cout << " - Rango h: " << min_h << " a " << max_h << "\n";
    std::cout << " - Rango mass: " << min_mass << " a " << max_mass << "\n";
    if (invalid_count > 0) {
        std::cout << " - Advertencia: " << invalid_count << " partículas con valores inválidos (h <= 0 o mass <= 0).\n";
    } else {
        std::cout << " - Todas las partículas tienen valores válidos.\n";
    }

    // Mostrar solo ejemplos
    std::cout << "\nEjemplos:\n";
    std::cout << " - Primera partícula -> ID: " << particles.front().id
              << " Pos(" << particles.front().pos[0] << ", " << particles.front().pos[1] << ")\n";
    std::cout << " - Última partícula -> ID: " << particles.back().id
              << " Pos(" << particles.back().pos[0] << ", " << particles.back().pos[1] << ")\n";
}
