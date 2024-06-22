#ifndef LEGENDRE_HPP
#define LEGENDRE_HPP

#include <vector>
#include <string>

#include "DataResult.h"

const std::string LegendreCosFileName = "python_legendre_coefficients_cos.txt";
const std::string LegendreSinFileName = "python_legendre_coefficients_sin.txt";

template<typename T>
std::vector<DataResult<T>> WorkLegendre(T x, std::function<T(T)> f, T result_x, bool isParallel, std::string fileCoefficients);

#include "../Legendre.cpp"

#endif // LEGENDRE_HPP
