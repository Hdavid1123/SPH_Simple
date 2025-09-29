#pragma once
#include <vector>
#include <array>      // necesario para std::array
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

// Construye la malla espacial
std::vector<Cell> buildCellGrid(double xmin, double xmax,
                                double ymin, double ymax,
                                double h, double kappa = 2.0);

// Calcula la caja mínima que contiene todas las partículas
void computeBoundingBox(const std::vector<Particle>& particles,
                        double& xmin, double& xmax,
                        double& ymin, double& ymax);

// 🚀 Función orquestadora: bounding box + malla + asignación + vecinos
std::vector<Cell> rebuildGridAndNeighbors(std::vector<Particle>& particles,
                                          double h, double kappa = 2.0);
