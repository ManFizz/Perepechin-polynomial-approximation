#ifndef TAYLOR_HPP
#define TAYLOR_HPP

#include <vector>
#include <string>

#include "DataResult.h"

template<typename T>
std::vector<DataResult<T>> WorkTaylor(T x, int maxCoefficient, std::function<T(T)> f);

#include "Taylor.cpp"

#endif //TAYLOR_HPP
