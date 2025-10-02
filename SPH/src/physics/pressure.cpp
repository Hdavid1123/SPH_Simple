#include "physics/pressure.h"
#include <cmath>

// EOS lineal: p = c^2 * rho
void computePressureEOS(std::vector<Particle>& particles, double c) {
    for (auto& pi : particles) {
        pi.soundVel = c;  // velocidad del sonido local (constante aquí)
        pi.pressure = c * c * pi.rho;
    }
}

// Tait EOS: p = c0^2 * rho0 * [ (rho / rho0)^gamma - 1 ]
void computePressureMonaghan(std::vector<Particle>& particles,
                         double rho0) {
    for (auto& pi : particles) {

        double c = 1480.0; // Velocidad del sonido en el agua  (m/s). Liu p. 400
        pi.soundVel = c;
        double gamma = 7.0;
        double b = 1.013e5; // Presión atmosférica en Pa
        
        pi.pressure = c * c * rho0 * (std::pow(pi.rho / rho0, gamma) - 1.0);
    }
}

void computePressureKorzani(std::vector<Particle>& particles,
                         double altura_columna, double g, double rho0) {
    
    for (auto& pi : particles) {
        // Velocidad del sonido en el agua  (m/s). Korzani pp. 3-4  
        double c_a = 10.0 * 2.0 * g * altura_columna; 
        pi.soundVel = c_a;
        //double c_a = 0.01;
        //pi.soundVel = c_a; // drho/rho
        double gamma = 7.0;
        double factor = rho0 * c_a * c_a / gamma;

        pi.pressure = factor * (std::pow((pi.rho / rho0), gamma) - 1.0);
    }
}
