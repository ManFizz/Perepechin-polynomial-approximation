#ifndef LEGENDRE_HPP
#define LEGENDRE_HPP

#include <vector>
#include <string>
#include <functional>

#include "DataResult.h"
#include "bignum.h"

template<typename T>
std::vector<DataResult<T>> WorkLegendre(T x, const int maxCoefficient, const int numPoints, std::function<T(T)> f, T result_x, bool isParallel);

#include "../Legendre.cpp"

#endif // LEGENDRE_HPP
