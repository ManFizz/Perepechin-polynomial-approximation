#include <vector>

#include "include/DataResult.h"
#include "include/helper.hpp"

template<typename T>
T computeFactorial(int n) {
    static std::vector<T> factorials = {1};
    if (n < static_cast<int>(factorials.size())) {
        return factorials[n];
    }
    T result = factorials.back();
    for (int i = factorials.size(); i <= n; ++i) {
        result *= i;
        factorials.push_back(result);
    }
    return result;
}

template<typename T>
T calculateCoefficientTaylor(int k, const std::function<T(T)>& f, T x0) {
    const T h = 1e-4;  // Adjusted for potential better stability
    std::vector<T> derivatives(k + 1, 0);
    derivatives[0] = f(x0);
    T factorial_i = 1;

    for (int i = 1; i <= k; ++i) {
        factorial_i *= i;  // compute factorial only once per loop
        T sum = 0;
        for (int j = 0; j <= i; ++j) {
            T factorial_j = computeFactorial<T>(j);
            T factorial_ij = computeFactorial<T>(i - j);
            T coef = ((j % 2 == 0) ? 1 : -1) * factorial_i / (factorial_j * factorial_ij);
            sum += coef * f(x0 + (i - 2 * j) * h);
        }
        derivatives[i] = sum / pow(2 * h, i);
    }

    return derivatives[k] / factorial_i;
}

template<typename T>
T approximateFunctionTaylor(const T& x, const std::vector<T>& coefficients, const T& x0) {
    T result = 0;
    T term = 1;
    for (size_t k = 0; k < coefficients.size(); ++k) {
        if (k > 0) {
            term *= (x - x0);
        }
        result += coefficients[k] * term;
    }
    return result;
}

template<typename T>
std::vector<DataResult<T>> WorkTaylor(const T& x, const int maxCoefficient, const std::function<T(T)>& f, const T& result_x) {
    T x0 = 0;
    std::vector<DataResult<T>> results;
    std::vector<T> coefficients;
    coefficients.reserve(maxCoefficient);

    for (int k = 0; k < maxCoefficient; ++k) {
        std::cout << "Taylor: " << k << std::endl;

        coefficients.push_back(calculateCoefficientTaylor(k, f, x0));

        auto start = std::chrono::high_resolution_clock::now();
        T approxValue = approximateFunctionTaylor(x, coefficients, x0);
        auto end = std::chrono::high_resolution_clock::now();
        DataResult<T>::AddData(results, approxValue, abs(result_x - approxValue), x, k, end - start);
    }

    return results;
}