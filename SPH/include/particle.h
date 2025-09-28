#pragma once
#include <array>
#include <vector>

const double RHO0 = 1000.0; // densidad de referencia (agua)

struct Particle {
    int id;
    std::array<double, 2> pos;     // posición (x,y)
    std::array<double, 2> vel;     // velocidad
    std::array<double, 2> accel;   // aceleración

    double mass;                   // masa
    double rho;                    // densidad actual
    double h;                      // longitud de suavizado
    double pressure;               // presión
    double soundVel;                      // velocidad del sonido local
    double internalE;                      // energía interna
    double dinternalE;                     // variación de energía

    std::vector<int> neighbors;    // índices de vecinos
    std::vector<double> dx;        // diferencia en x con vecinos
    std::vector<double> dy;        // diferencia en y con vecinos
    std::vector<double> r;         // distancia con vecinos
    std::vector<double> W;         // valores del kernel
    std::vector<double> dWx;       // derivadas parciales
    std::vector<double> dWy;

    int type;                      // tipo de partícula (fluido, frontera, etc.)
};
