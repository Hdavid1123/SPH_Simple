#pragma once
#include <vector>
#include "core/particle.h"

// EOS lineal: p = c^2 * rho
void computePressureEOS(std::vector<Particle>& particles, double c);

// Tait EOS: p = c0^2 * rho0 * [ (rho / rho0)^gamma - 1 ]
void computePressureTait(std::vector<Particle>& particles,
                         double rho0, double c0, double gamma);
