#pragma once
#include <vector>
#include "core/particle.h"

// EOS lineal: p = c^2 * rho
void computePressureEOS(std::vector<Particle>& particles, double c);

// Tait (Liu 2003) EOS: p = c0^2 * rho0 * [ (rho / rho0)^gamma - 1 ]
void computePressureMonaghan(std::vector<Particle>& particles,
                         double rho0);

// Tait (Korzani 2014)
void computePressureKorzani(std::vector<Particle>& particles,
                         double altura_columna, double g, double rho0);

