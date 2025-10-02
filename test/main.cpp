#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <vector>
#include <unordered_set>
#include <cmath>

// =======================
// Definiciones mínimas
// =======================
struct Particle {
    int id;
    double h;
    int type;
    double mass;
    std::array<double,2> pos;
    std::array<double,2> vel;
    std::array<double,2> accel;
    double rho;
    double pressure;
    double soundVel;
    double internalE;
    double dinternalE;

    // vecinos
    std::vector<int> neighbors;
    std::vector<double> dx, dy, r, W, dWx, dWy;
};

struct Cell {
    int id;
    std::vector<double> center;
    std::vector<int> particles;
    std::vector<int> neighborCells;
};

// =======================
// Lectura de partículas
// =======================
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
        if (line.empty()) continue;
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

    std::cout << "Resumen de validación:\n";
    std::cout << " - Partículas leídas: " << particles.size() << "\n";

    std::cout << "\nEjemplos:\n";
    std::cout << " - Primera partícula -> ID: " << particles.front().id
              << " Pos(" << particles.front().pos[0] << ", " << particles.front().pos[1] << ")\n";
    std::cout << " - Última partícula -> ID: " << particles.back().id
              << " Pos(" << particles.back().pos[0] << ", " << particles.back().pos[1] << ")\n";
}

// =======================
// Funciones de vecinos
// =======================
void freeParticleNeighbors(std::vector<Particle>& particles) {
    for (auto& p : particles) {
        std::vector<int>().swap(p.neighbors);
        std::vector<double>().swap(p.dx);
        std::vector<double>().swap(p.dy);
        std::vector<double>().swap(p.r);
        std::vector<double>().swap(p.W);
        std::vector<double>().swap(p.dWx);
        std::vector<double>().swap(p.dWy);
    }
}

void assignParticlesToCells(std::vector<Cell>& cells,
                            std::vector<Particle>& particles,
                            double h, double kappa) {
    double cell_size = kappa * h;

    for (auto& cell : cells) {
        cell.particles.clear();
        for (size_t j = 0; j < particles.size(); j++) {
            double dx = std::abs(particles[j].pos[0] - cell.center[0]);
            double dy = std::abs(particles[j].pos[1] - cell.center[1]);

            if (dx < 0.5 * cell_size + 1e-5 && dy < 0.5 * cell_size + 1e-5) {
                cell.particles.push_back(j);
            }
        }
    }
}

void findNeighbors(std::vector<Cell>& cells,
                   std::vector<Particle>& particles,
                   double kappa) {
    for (auto& cell : cells) {
        for (int pid : cell.particles) {
            Particle& pi = particles[pid];

            pi.neighbors.clear();
            pi.dx.clear(); pi.dy.clear(); pi.r.clear();
            pi.W.clear(); pi.dWx.clear(); pi.dWy.clear();

            std::unordered_set<int> uniqueNeighbors;

            for (int neighborCellId : cell.neighborCells) {
                Cell& neighCell = cells[neighborCellId];
                for (int pj_id : neighCell.particles) {
                    if (pj_id == pid) continue;

                    Particle& pj = particles[pj_id];
                    double dx = pi.pos[0] - pj.pos[0];
                    double dy = pi.pos[1] - pj.pos[1];
                    double r = std::sqrt(dx*dx + dy*dy);
                    double h_ij = 0.5 * (pi.h + pj.h);

                    if (r < kappa * h_ij) {
                        if (uniqueNeighbors.insert(pj_id).second) {
                            pi.neighbors.push_back(pj_id);
                            pi.dx.push_back(dx);
                            pi.dy.push_back(dy);
                            pi.r.push_back(r);
                            pi.W.push_back(1.0 / (1.0 + r)); // placeholder
                            pi.dWx.push_back(-dx);
                            pi.dWy.push_back(-dy);
                        }
                    }
                }
            }
        }
    }
}

std::vector<Cell> buildCellGrid(double xmin, double xmax,
                                double ymin, double ymax,
                                double h, double kappa) {
    double cell_size = kappa * h;
    int nx = static_cast<int>(std::ceil((xmax - xmin) / cell_size));
    int ny = static_cast<int>(std::ceil((ymax - ymin) / cell_size));

    std::vector<Cell> cells;
    cells.reserve(nx * ny);

    auto cellIndex = [nx](int i, int j) { return j * nx + i; };

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            Cell cell;
            cell.id = cellIndex(i, j);
            cell.center = { xmin + (i + 0.5) * cell_size,
                            ymin + (j + 0.5) * cell_size };

            for (int dj = -1; dj <= 1; dj++) {
                for (int di = -1; di <= 1; di++) {
                    int ni = i + di;
                    int nj = j + dj;
                    if (ni >= 0 && ni < nx && nj >= 0 && nj < ny) {
                        cell.neighborCells.push_back(cellIndex(ni, nj));
                    }
                }
            }
            cells.push_back(std::move(cell));
        }
    }
    return cells;
}

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

std::vector<Cell> rebuildGridAndNeighbors(std::vector<Particle>& particles,
                                          double h, double kappa) {
    double xmin, xmax, ymin, ymax;
    computeBoundingBox(particles, xmin, xmax, ymin, ymax);

    auto cells = buildCellGrid(xmin - h, xmax + h,
                               ymin - h, ymax + h,
                               h, kappa);

    assignParticlesToCells(cells, particles, h, kappa);
    findNeighbors(cells, particles, kappa);

    return cells;
}

// =======================
// MAIN de prueba
// =======================
int main(int argc, char* argv[]) {
    std::string filename = (argc > 1) ? argv[1] : "particulas.txt";

    auto particles = readParticlesFromFile(filename);
    validateParticles(particles);

    double h = 0.1, kappa = 2.0;
    auto cells = rebuildGridAndNeighbors(particles, h, kappa);

    std::cout << "\nVecinos iniciales:\n";
    for (auto& p : particles) {
        std::cout << "Particle " << p.id
                  << " tiene " << p.neighbors.size() << " vecinos\n";
    }

    freeParticleNeighbors(particles);

    std::cout << "\nTras liberar memoria:\n";
    for (auto& p : particles) {
        std::cout << "Particle " << p.id
                  << " vecinos=" << p.neighbors.size()
                  << " dx=" << p.dx.size()
                  << " W=" << p.W.size() << "\n";
    }

    return 0;
}
