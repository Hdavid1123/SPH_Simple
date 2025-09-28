#pragma once
#include <vector>
#include "core/particle.h"

// Avanza posiciones y energía interna medio paso
void drift(std::vector<Particle>& particles, double dt);

// Avanza velocidades un paso completo con aceleraciones
void kick(std::vector<Particle>& particles, double dt);
