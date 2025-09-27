#pragma once
#include <vector>
#include "particle.h"

// Calcula el bounding box dinámico del conjunto de partículas
void computeBoundingBox(const std::vector<Particle>& particles,
                        double& xmin, double& xmax,
                        double& ymin, double& ymax);

// Test de vecinos: escribe en NN_test.output
void test_NN(const std::vector<Particle>& particles, int nTests);
