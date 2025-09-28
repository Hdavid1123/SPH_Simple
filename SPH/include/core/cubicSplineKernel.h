#pragma once
#include <array>
#include <stdexcept>
#include <cmath>


enum Dimension { ONE_D = 1, TWO_D = 2, THREE_D = 3 };

double cubicSplineKernel(double r, double h, Dimension dim);

std::array<double,2> dCubicSplineKernel(double r, double dx, double dy, double h);

