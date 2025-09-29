#pragma once
#include <vector>
#include <string>
#include "core/particle.h"
#include "core/linkedList.h"

// Inicializa una fase de simulación:
// - reconstruye malla y vecinos
// - calcula densidad y presión
// - imprime datos y guarda test
std::vector<Cell> initializePhase(
    std::vector<Particle>& particles, double h_ref, double kappa, double c,
    const std::string& label, const std::string& testFile
);


// Test de vecinos: escribe en NN_test.output
void test_NN(const std::vector<Particle>& particles, int nTests,
             const std::string& filename);