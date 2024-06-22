#ifndef CHEBYSHEV_HPP
#define CHEBYSHEV_HPP

#include <vector>
#include <string>

#include "DataResult.h"

const std::string ChebyshevCosFileName = "python_chebyshev_coefficients_cos.txt";
const std::string ChebyshevSinFileName = "python_chebyshev_coefficients_sin.txt";

template<typename T>
std::vector<DataResult<T>> WorkChebyshev(T x, std::function<T(T)> f, T result_x, bool isParallel, std::string fileCoefficients);

#include "../Chebyshev.cpp"

#endif // CHEBYSHEV_HPP
