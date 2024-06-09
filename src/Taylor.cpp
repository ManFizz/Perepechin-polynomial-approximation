#include <vector>

#include "include/DataResult.h"
#include "include/helper.hpp"

bigfloat_t computeFactorial(int n) {
    static std::vector<bigfloat_t> factorials = {1};
    if (n < static_cast<int>(factorials.size())) {
        return factorials[n];
    }
    bigfloat_t result = factorials.back();
    for (int i = factorials.size(); i <= n; ++i) {
        result *= i;
        factorials.push_back(result);
    }
    return result;
}

bigfloat_t calculateCoefficientTaylor(int k, const std::function<bigfloat_t(bigfloat_t)>& f, bigfloat_t& x0) {
    const bigfloat_t h = 1e-4;  // Adjusted for potential better stability
    std::vector<bigfloat_t> derivatives(k + 1, 0);
    derivatives[0] = f(x0);
    bigfloat_t factorial_i = 1;

    for (int i = 1; i <= k; ++i) {
        factorial_i *= i;  // compute factorial only once per loop
        bigfloat_t sum = 0;
        for (int j = 0; j <= i; ++j) {
            bigfloat_t factorial_j = computeFactorial(j);
            bigfloat_t factorial_ij = computeFactorial(i - j);
            bigfloat_t coef = ((j % 2 == 0) ? 1 : -1) * factorial_i / (factorial_j * factorial_ij);
            sum += coef * f(x0 + (i - 2 * j) * h);
        }
        derivatives[i] = sum / pow(2 * h, i);
    }

    return derivatives[k] / factorial_i;
}

bigfloat_t approximateFunctionTaylor(const bigfloat_t& x, const std::vector<bigfloat_t>& coefficients, const bigfloat_t& x0) {
    bigfloat_t result = 0;
    bigfloat_t term = 1;
    for (size_t k = 0; k < coefficients.size(); ++k) {
        if (k > 0) {
            term *= (x - x0);
        }
        result += coefficients[k] * term;
    }
    return result;
}

std::vector<DataResult<bigfloat_t>> WorkTaylor(bigfloat_t& x, int maxCoefficient, std::function<bigfloat_t(bigfloat_t)> f, bigfloat_t& result_x) {
    bigfloat_t x0 = 0;
    std::vector<DataResult<bigfloat_t>> results;
    std::vector<bigfloat_t> coefficients;
    coefficients.reserve(maxCoefficient);

    for (int k = 0; k < maxCoefficient; ++k) {
        coefficients.push_back(calculateCoefficientTaylor(k, f, x0));

        auto start = std::chrono::high_resolution_clock::now();
        bigfloat_t approxValue = approximateFunctionTaylor(x, coefficients, x0);
        auto end = std::chrono::high_resolution_clock::now();
        DataResult<bigfloat_t>::AddData(results, approxValue, abs(result_x - approxValue), x, k, end - start);
    }

    return results;
}