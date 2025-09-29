#include "utils/output.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <stdexcept>

void printState(const std::vector<Particle>& particles, int step) {
    namespace fs = std::filesystem;
    fs::create_directory("output");  // crea carpeta si no existe

    std::ostringstream filename;
    filename << "output/state_" << std::setw(4) << std::setfill('0') << step << ".txt";

    std::ofstream out(filename.str());
    if (!out) {
        throw std::runtime_error("No se pudo abrir archivo " + filename.str());
    }

    for (const auto& p : particles) {
        out << p.id << " "
            << p.pos[0] << " " << p.pos[1] << " "
            << p.vel[0] << " " << p.vel[1] << " "
            << p.accel[0] << " " << p.accel[1] << " "
            << p.rho << " "
            << p.mass << " "
            << p.pressure << " "
            << p.h << " "
            << p.type
            << "\n";
    }
}
