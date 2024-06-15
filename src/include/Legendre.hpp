#ifndef LEGENDRE_HPP
#define LEGENDRE_HPP

#include <vector>
#include <string>

#include "DataResult.h"

const std::string LegendreCosFileName = "legendre_coefficients_cos.txt";
const std::string LegendreSinFileName = "legendre_coefficients_sin.txt";

template<typename T>
std::vector<DataResult<T>> WorkLegendre(T x, int maxCoefficient, int numPoints, std::function<T(T)> f, T result_x, bool isParallel, std::string fileCoefficients);

#include "../Legendre.cpp"

#endif // LEGENDRE_HPP
