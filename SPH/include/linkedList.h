#pragma once
#include <vector>
#include "particle.h"

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
                            double h);

// Encuentra vecinos usando la lista enlazada de celdas
void findNeighbors(std::vector<Cell>& cells,
                   std::vector<Particle>& particles,
                   double kappa);

// Construye la malla espacial y asigna partículas a las celdas
std::vector<Cell> buildCellGrid(double xmin, double xmax,
                                double ymin, double ymax,
                                double h);