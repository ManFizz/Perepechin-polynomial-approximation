#include <cmath>
#include <vector>
#include <functional>

#include "DataResult.h"

template<typename T>
T chebyshevPolynomial(int n, T x) {
    if (n == 0) return T(1);
    if (n == 1) return x;

    return 2 * x * chebyshevPolynomial(n - 1, x) - chebyshevPolynomial(n - 2, x);
}

template<typename T>
T calculateCoefficientChebyshev(int k, int numPoints, std::function<T(T)> f) {
    T coefficient = T(0);
    T weight = (k == 0) ? T(1.0 / numPoints) : T(2.0 / numPoints);

    for (int n = 0; n < numPoints; ++n) {
        T x = std::cos(M_PI * (T(n) + T(0.5)) / T(numPoints));
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
std::vector<DataResult<T>> WorkChebyshev(const T x, const int maxCoefficient, const int numPoints, std::function<T(T)> f) {
    std::vector<DataResult<T>> results;
    std::vector<T> coefficients = {};
    for(int k = 0; k < maxCoefficient; k++) {
        T result = calculateCoefficientChebyshev<T>(k, numPoints, f);
        coefficients.push_back(result);

        T reslut = approximateFunctionChebyshev<T>(x, coefficients);
        DataResult<T>::AddData(results, reslut, x, k);
    }
    return results;
}