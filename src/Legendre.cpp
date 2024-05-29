#include <vector>
#include <cmath>
#include <functional>

#include "include/DataResult.h"

template<typename T>
T legendrePolynomial(int n, T x) {
    if (n == 0) return T(1);
    if (n == 1) return x;
    return ((2 * n - 1) * x * legendrePolynomial(n - 1, x) - (n - 1) * legendrePolynomial(n - 2, x)) / T(n);
}

template<typename T>
T trapezoidalRule(std::function<T(T)> func, T a, T b, int n) {
    T h = (b - a) / n;
    T sum = func(a) / 2.0 + func(b) / 2.0;
    for (int i = 1; i < n; i++) {
        T x_i = a + i * h;
        sum += func(x_i);
    }
    return sum * h;
}

template<typename T>
T integrandFunction(T x, int k, std::function<T(T)> f) {
    return f(x) * legendrePolynomial(k, x);
}

template<typename T>
T calculateCoefficientLegendre(int k, int numPoints, std::function<T(T)> f) {
    T a = -1.0, b = 1.0;
    std::function<T(T)> func = [k, f](T x) -> T {
        return integrandFunction(x, k, f);
    };
    T coefficient = (2 * k + 1) / 2.0 * trapezoidalRule(func, a, b, numPoints);
    return coefficient;
}


template<typename T>
T approximateFunctionLegendre(T x, const std::vector<T>& coefficients) {
    T sum = T(0.0);
    for (size_t k = 0; k < coefficients.size(); ++k) {
        sum += coefficients[k] * legendrePolynomial(static_cast<int>(k), x);
    }
    return sum;
}

template<typename T>
std::vector<DataResult<T>> WorkLegendre(const T x, const int maxCoefficient, const int numPoints, std::function<T(T)> f, const T result_x) {
    std::vector<DataResult<T>> results;
    std::vector<T> coefficients = {};
    for (int k = 0; k < maxCoefficient; k++) {
        T coefficient = calculateCoefficientLegendre(k, numPoints, f);
        coefficients.push_back(coefficient);

        T result = approximateFunctionLegendre(x, coefficients);
        DataResult<T>::AddData(results, std::abs(result_x - result), x, k);
    }
    return results;
}