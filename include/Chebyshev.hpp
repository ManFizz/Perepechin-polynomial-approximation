#ifndef CHEBYSHEV_HPP
#define CHEBYSHEV_HPP

#include <vector>
#include <string>

#include "DataResult.h"

template<typename T>
std::vector<DataResult<T>> WorkChebyshev(T x, int maxCoefficient, int numPoints, std::function<T(T)> f);

#include "Chebyshev.cpp"

#endif // CHEBYSHEV_HPP
