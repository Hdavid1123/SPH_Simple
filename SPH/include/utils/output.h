#ifndef OUTPUT_H
#define OUTPUT_H

#include <vector>
#include "../core/particle.h"  // Asegúrate que aquí está definida la struct/class Particle

// Guarda el estado de todas las partículas en un archivo ./output/state_XXXX.txt
void printState(const std::vector<Particle>& particles, int step);

#endif
