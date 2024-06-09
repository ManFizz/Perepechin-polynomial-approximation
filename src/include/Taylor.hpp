#ifndef TAYLOR_HPP
#define TAYLOR_HPP

#include <vector>
#include <string>

#include "DataResult.h"
#include "bignum.h"

std::vector<DataResult<bigfloat_t>> WorkTaylor(bigfloat_t& x, int maxCoefficient, std::function<bigfloat_t(bigfloat_t)> f, bigfloat_t& result_x);

#include "../Taylor.cpp"

#endif //TAYLOR_HPP
