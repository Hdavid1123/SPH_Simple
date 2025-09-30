#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include "core/cubicSplineKernel.h"

double cubicSplineKernel(double r, double h, Dimension dim) {
    if (r < 0 || h <= 0) throw std::runtime_error("r>=0 and h>0 required");

    double q = r / h;
    double alpha;

    switch(dim) {
        case ONE_D:   alpha = 2.0 / 3.0; break;
        case TWO_D:   alpha = 10.0 / (7.0 * M_PI); break;
        case THREE_D: alpha = 1.0 / M_PI; break;
        default:      throw std::runtime_error("Dimensión no soportada");
    }

    alpha /= std::pow(h, static_cast<int>(dim));

    double w = 0.0;
    if (q >= 0.0 && q < 1.0) {
        w = 1.0 - 1.5*q*q + 0.75*q*q*q;
    } else if (q >= 1.0 && q < 2.0) {
        w = 0.25 * std::pow(2.0 - q, 3);
    } else if (q > 2.0){
        w = 0.0;
    }
    return alpha * w;
}

// Derivadas direccionales dWx y dWy (solo en 2D por ahora)
std::array<double,2> dCubicSplineKernel(double r, double dx, double dy, double h) {
    if (r <= 0 || h <= 0) return {0.0, 0.0}; // evita división por cero

    double R = r / h;
    double alpha = 15.0 / (7.0 * M_PI * h * h); // normalización 2D
    double factor = 0.0;

    if (R >= 0.0 && R < 1.0) {
        factor = alpha * (-2.0 + 1.5*R) / (h*h);
    } else if (R>= 1.0 && R < 2.0) {
        factor = alpha * (-0.5*(2.0-R)*(2.0-R)) / (h*h*R);
    } else {
        factor = 0.0;
    }

    double dWx = factor * dx;
    double dWy = factor * dy;

    return {dWx, dWy};
}

// Función de prueba del kernel cúbico
void testKernel(double h, Dimension dim) {
    std::ofstream fKernelTest("test_results/kernel_test.output");
    if (!fKernelTest.is_open()) {
        throw std::runtime_error("No se pudo abrir kernel_test.output");
    }

    fKernelTest << std::setw(16) << "r"
                << std::setw(16) << "W"
                << std::setw(16) << "dWx"
                << std::setw(16) << "dWy" << "\n";

    for (double r = -3.0; r < 3.0; r += 0.1) {
        double rr = std::fabs(r);

        // Evaluamos el kernel en 1D, 2D o 3D según corresponda
        double w = cubicSplineKernel(rr, h, dim);

        // En este ejemplo tomamos dirección en x para derivada
        std::array<double,2> dw = dCubicSplineKernel(rr, r, 0.0, h);

        fKernelTest << std::setw(16) << std::fixed << std::setprecision(10) << r
                    << std::setw(16) << w
                    << std::setw(16) << dw[0]
                    << std::setw(16) << dw[1]
                    << "\n";
            
    }

    std::cout << "Prueba del kernel completada\n.";
    std::cout << "Resultados en test_results/kernel_test.output\n";

    fKernelTest.close();
}
