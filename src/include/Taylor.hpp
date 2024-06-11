#ifndef TAYLOR_HPP
#define TAYLOR_HPP

#include <vector>
#include <string>

#include "DataResult.h"
#include "bignum.h"

std::vector<DataResult<bigfloat_t>> WorkTaylor(bigfloat_t& x, int maxCoefficient, bigfloat_t& result_x, std::function<bigfloat_t(bigfloat_t&, int)> f);

bigfloat_t approximateCosTaylor(bigfloat_t& x, int terms);

bigfloat_t approximateSinTaylor(bigfloat_t& x, int terms);

#include "../Taylor.cpp"

#endif //TAYLOR_HPP
