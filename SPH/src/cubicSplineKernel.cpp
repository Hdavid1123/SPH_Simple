#include "cubicSplineKernel.h"

double cubicSplineKernel(double r, double h, Dimension dim) {
    if (r < 0 || h <= 0) throw std::runtime_error("r>=0 and h>0 required");

    double q = r / h;
    double sigma;

    switch(dim) {
        case ONE_D:   sigma = 2.0 / 3.0; break;
        case TWO_D:   sigma = 10.0 / (7.0 * M_PI); break;
        case THREE_D: sigma = 1.0 / M_PI; break;
        default:      throw std::runtime_error("Dimensión no soportada");
    }

    sigma /= std::pow(h, static_cast<int>(dim));

    double w = 0.0;
    if (q >= 0.0 && q < 1.0) {
        w = 1.0 - 1.5*q*q + 0.75*q*q*q;
    } else if (q >= 1.0 && q < 2.0) {
        w = 0.25 * std::pow(2.0 - q, 3);
    }

    return sigma * w;
}

// Derivadas direccionales dWx y dWy (solo en 2D por ahora)
std::array<double,2> dCubicSplineKernel(double r, double dx, double dy, double h) {
    if (r <= 0 || h <= 0) return {0.0, 0.0}; // evita división por cero

    double R = r / h;
    double alpha = 15.0 / (7.0 * M_PI * h * h); // normalización 2D
    double factor = 0.0;

    if (R < 1.0) {
        factor = alpha * (-2.0 + 1.5*R) / (h*h);
    } else if (R < 2.0) {
        factor = alpha * (-0.5*(2.0-R)*(2.0-R)) / (h*h*R);
    } else {
        return {0.0, 0.0};
    }

    double dWx = factor * dx;
    double dWy = factor * dy;

    return {dWx, dWy};
}
