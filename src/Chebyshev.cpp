#include <cmath>
#include <vector>

#include "include/DataResult.h"
#include "include/helper.hpp"

template<typename T>
T chebyshevPolynomial(int n, T x) {
    if (n == 0) return T(1);
    if (n == 1) return x;

    T Tn_2 = T(1);
    T Tn_1 = x;
    T Tn;

    for (int i = 2; i <= n; ++i) {
        Tn = 2 * x * Tn_1 - Tn_2;
        Tn_2 = Tn_1;
        Tn_1 = Tn;
    }

    return Tn;
}

template<typename T>
T calculateCoefficientChebyshev(int k, int numPoints, std::function<T(T)> f) {
    T coefficient = T(0);
    T weight = (k == 0) ? T(1.0 / numPoints) : T(2.0 / numPoints);

    for (int n = 0; n < numPoints; ++n) {
        T x = cos(bigfloat_t::pi * (T(n) + T(0.5)) / T(numPoints), 60);
        coefficient += weight * f(x) * chebyshevPolynomial<T>(k, x);
    }

    return coefficient;
}

template<typename T>
T approximateFunctionChebyshev(T x, const std::vector<T>& coefficients) {
    T sum = T(0.0);
    for (size_t k = 0; k < coefficients.size(); ++k) {
        sum += coefficients[k] * chebyshevPolynomial<T>((int)k, x);
    }
    return sum;
}

template<typename T>
std::vector<DataResult<T>> WorkChebyshev(T x, int maxCoefficient, int numPoints, std::function<T(T)> f, T result_x) {
    std::vector<DataResult<T>> results;
    std::vector<T> coefficients = {};

    std::cout << "Chebyshev:" << std::endl;
    for(int k = 0; k < maxCoefficient; k++) {

        coefficients.push_back(calculateCoefficientChebyshev(k, numPoints, f));

        auto start = std::chrono::high_resolution_clock::now();
        T approxValue = approximateFunctionChebyshev<T>(x, coefficients);
        auto end = std::chrono::high_resolution_clock::now();

        auto result = DataResult<T>(approxValue, abs(result_x - approxValue), x, k, end - start);
        results.emplace_back(result);

        std::cout << result << std::endl;
    }
    return results;
}