#include <unordered_set>
#include <cmath>
#include <limits>  // para std::numeric_limits
#include "core/cubicSplineKernel.h"
#include "core/linkedList.h"

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

                            double Wval = cubicSplineKernel(r, h_ij, TWO_D);
                            auto dWval = dCubicSplineKernel(r, dx, dy, h_ij);

                            pi.W.push_back(Wval);
                            pi.dWx.push_back(dWval[0]);
                            pi.dWy.push_back(dWval[1]);
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

    auto cellIndex = [nx](int i, int j) {
        return j * nx + i;  // row-major
    };

    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            Cell cell;
            cell.id = cellIndex(i, j);
            cell.center = { xmin + (i + 0.5) * cell_size,
                            ymin + (j + 0.5) * cell_size };

            // vecinos: 8 alrededor + la propia
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
