#pragma once
#include <vector>
#include "core/particle.h"

// Calcula aceleraciones por ecuaciones de Navier-Stokes
void computeNavierStokes(std::vector<Particle>& particles);
