#pragma once
#include <array>
#include <stdexcept>
#include <cmath>


enum Dimension { ONE_D = 1, TWO_D = 2, THREE_D = 3 };

// Función de kernel cúbico spline
double cubicSplineKernel(double r, double h, Dimension dim);

// Derivada del kernel cúbico spline
std::array<double,2> dCubicSplineKernel(double r, double dx, double dy, double h);

// Función de prueba para verificar la implementación del kernel y su derivada
void testKernel(double h, Dimension dim);


