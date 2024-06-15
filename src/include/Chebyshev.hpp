#ifndef CHEBYSHEV_HPP
#define CHEBYSHEV_HPP

#include <vector>
#include <string>

#include "DataResult.h"

const std::string ChebyshevCosFileName = "chebyshev_coefficients_cos.txt";
const std::string ChebyshevSinFileName = "chebyshev_coefficients_sin.txt";

template<typename T>
std::vector<DataResult<T>> WorkChebyshev(T x, int maxCoefficient, int numPoints, std::function<T(T)> f, T result_x, bool isParallel, std::string fileCoefficients);

#include "../Chebyshev.cpp"

#endif // CHEBYSHEV_HPP
