#pragma once
#include <vector>
#include "core/particle.h"

// Estructura para celdas de la malla espacial
struct Cell {
    int id;
    std::array<double,2> center;
    std::vector<int> particles;
    std::vector<int> neighborCells;
};

// Asigna partículas a las celdas
void assignParticlesToCells(std::vector<Cell>& cells,
                            std::vector<Particle>& particles,
                            double h, double kappa = 2.0);

// Encuentra vecinos usando la lista enlazada de celdas
void findNeighbors(std::vector<Cell>& cells,
                   std::vector<Particle>& particles,
                   double kappa = 2.0);

// Construye la malla espacial y asigna partículas a las celdas
std::vector<Cell> buildCellGrid(double xmin, double xmax,
                                double ymin, double ymax,
                                double h, double kappa = 2.0);