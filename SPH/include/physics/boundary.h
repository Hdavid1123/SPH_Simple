#pragma once
#include <vector>
#include "core/particle.h"

// Interacción con las partículas de frontera
void boundaryInteraction(std::vector<Particle>& particles);

// Versión con fuerza cuadrática
void boundaryInteractionCuadratic(std::vector<Particle>& particles);
