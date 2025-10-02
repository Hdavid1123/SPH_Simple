#pragma once
#include <vector>
#include "core/particle.h"

// Calcula todas las aceleraciones relevantes (viscosidad, Navier-Stokes, frontera, gravedad)
void acceleration(std::vector<Particle>& particles, double g, double dx);

// Viscosidad artificial estilo Monaghan-Gingold
void viscosity(std::vector<Particle>& particles);

// Interacción con frontera
void boundaryInteraction(std::vector<Particle>& particles, double dx);

// Aceleración gravitacional
void applyGravity(std::vector<Particle>& particles, double g);

// Compute Navier-Stokes (llama viscosidad y gravedad)
void computeNavierStokes(std::vector<Particle>& particles);
