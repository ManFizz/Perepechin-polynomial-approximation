#include <vector>
#include <iostream>
#include <omp.h>

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

template<typename T>
T approximateFunctionLegendre(const T& x, const std::vector<T>& coefficients) {
    T sum = T(0.0);
    for (size_t k = 0; k < coefficients.size(); ++k) {
        sum += coefficients[k] * legendrePolynomial((int)(k), x);
    }
    return sum;
}

template<typename T>
T approximateFunctionLegendreOMP(const T& x, const std::vector<T>& coefficients) {
    int num_threads = omp_get_max_threads();
    std::vector<T> local_sums(num_threads, T(0.0));

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for
        for (size_t k = 0; k < coefficients.size(); ++k) {
            local_sums[thread_id] += coefficients[k] * legendrePolynomial((int)k, x);
        }
    }

    T sum = T(0.0);
    for (T local_sum : local_sums) {
        sum += local_sum;
    }

    return sum;
}

template<typename T>
std::vector<DataResult<T>> WorkLegendre(T x, const int maxCoefficient, const int numPoints, std::function<T(T)> f, T result_x, bool isParallel) {
    std::vector<DataResult<T>> results;
    std::vector<T> xValues(numPoints);
    std::vector<T> yValues(numPoints);
    T a = -1.0, b = 1.0;
    T h = (b - a) / (numPoints - 1);

    for (int i = 0; i < numPoints; ++i) {
        xValues[i] = a + i * h;
        yValues[i] = f(xValues[i]);
    }

    std::cout << "Legendre:"<< std::endl;
    for (int k = 1; k <= maxCoefficient; ++k) {

        std::vector<T> coefficients = fitLeastSquares(xValues, yValues, k, legendrePolynomial<T>);
        T approxValue;
        auto start = std::chrono::high_resolution_clock::now();
        if (isParallel) {
            approxValue = approximateFunctionLegendreOMP(x, coefficients);
        } else {
            approxValue = approximateFunctionLegendre(x, coefficients);
        }
        auto end = std::chrono::high_resolution_clock::now();

        auto result = DataResult<T>(approxValue, abs(result_x - approxValue), x, k, end - start);
        results.emplace_back(result);

        std::cout << result << std::endl;
    }

    return results;
}