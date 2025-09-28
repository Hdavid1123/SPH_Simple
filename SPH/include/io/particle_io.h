#pragma once
#include <vector>
#include <string>
#include "core/particle.h"

// Lee archivo de condiciones iniciales y devuelve un vector de partículas
std::vector<Particle> readParticlesFromFile(const std::string& filename);

// Valida y resume los datos de las partículas leídas
void validateParticles(const std::vector<Particle>& particles);
