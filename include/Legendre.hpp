#ifndef LEGENDRE_HPP
#define LEGENDRE_HPP

#include <vector>
#include <string>
#include <functional>

#include "DataResult.h"

template<typename T>
std::vector<DataResult<T>> WorkLegendre(T x, int maxCoefficient, int numPoints, std::function<T(T)> f, T result_x);

#include "Legendre.cpp"

#endif // LEGENDRE_HPP
