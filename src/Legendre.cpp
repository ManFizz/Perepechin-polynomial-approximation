#include <vector>
#include <functional>
#include <iostream>

#include "include/DataResult.h"
#include "include/helper.hpp"

template<typename T>
T legendrePolynomial(int n, T x) {
    if (n == 0) return T(1);
    if (n == 1) return x;

    T P0 = T(1);
    T P1 = x;
    T Pn = T(0);

    for (int i = 2; i <= n; ++i) {
        Pn = ((2 * i - 1) * x * P1 - (i - 1) * P0) / T(i);
        P0 = P1;
        P1 = Pn;
    }

    return Pn;
}

// Аппроксимация функции Лежандра
template<typename T>
T approximateFunctionLegendre(const T& x, const std::vector<T>& coefficients) {
    T sum = T(0.0);
    for (size_t k = 0; k < coefficients.size(); ++k) {
        sum += coefficients[k] * legendrePolynomial(static_cast<int>(k), x);
    }
    return sum;
}

template<typename T>
std::vector<DataResult<T>> WorkLegendre(T x, const int maxCoefficient, const int numPoints, std::function<T(T)> f, T result_x) {
    std::vector<DataResult<T>> results;
    std::vector<T> xValues(numPoints);
    std::vector<T> yValues(numPoints);
    T a = -1.0, b = 1.0;
    T h = (b - a) / (numPoints - 1);

    for (int i = 0; i < numPoints; ++i) {
        xValues[i] = a + i * h;
        yValues[i] = f(xValues[i]);
    }

    for (int k = 1; k <= maxCoefficient; ++k) {
        std::vector<T> coefficients = fitLeastSquares(xValues, yValues, k, legendrePolynomial<T>);
        T result = approximateFunctionLegendre(x, coefficients);
        DataResult<T>::AddData(results, abs(result_x - result), x, k);
    }

    return results;
}