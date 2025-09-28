#pragma once
#include <vector>
#include "particle.h"

// Calcula densidad SPH estándar: rho_i = sum_j m_j W(r_ij, h_ij)
void computeDensity(std::vector<Particle>& particles);

// Actualiza h_i adaptativo según rho_i y rho0
void updateSmoothingLength(std::vector<Particle>& particles, double rho0);
