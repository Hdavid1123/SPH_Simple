#include "physics/integrator.h"

void drift(std::vector<Particle>& particles, double dt) {
    for (auto& p : particles) {
        if (p.type != 0) continue; // solo fluido
        p.pos[0] += 0.5 * dt * p.vel[0];
        p.pos[1] += 0.5 * dt * p.vel[1];
        p.internalE += 0.5 * dt * p.dinternalE;
    }
}

void kick(std::vector<Particle>& particles, double dt) {
    for (auto& p : particles) {
        if (p.type != 0) continue; // solo fluido
        p.vel[0] += dt * p.accel[0];
        p.vel[1] += dt * p.accel[1];
    }
}
