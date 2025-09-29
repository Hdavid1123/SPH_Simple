#pragma once
#include <vector>
#include "core/particle.h"

// Muestra un resumen simple de densidad y presión
void printDensityPressure(const std::vector<Particle>& particles, int nPrint = 5);

// Muestra un resumen simple de aceleraciones
void printAccelerations(const std::vector<Particle>& particles, int nPrint = 5);

// Imprime información de verificación sobre la primera partícula fluida encontrada
void verifyFirstFluidParticle(const std::vector<Particle>& particles);
