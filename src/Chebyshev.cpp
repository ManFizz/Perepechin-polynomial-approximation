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
T calculateCoefficientChebyshevOMP(int k, int numPoints, std::function<T(T)> f) {
    std::vector<T> local_sums(omp_get_max_threads(), T(0));

    T weight = (k == 0) ? T(1.0 / numPoints) : T(2.0 / numPoints);

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for
        for (int n = 0; n < numPoints; ++n) {
            T x = cos(bigfloat_t::pi * (T(n) + T(0.5)) / T(numPoints), 60);
            local_sums[thread_id] += weight * f(x) * chebyshevPolynomial<T>(k, x);
        }
    }

    T coefficient = T(0);
    for (auto sum : local_sums) {
        coefficient += sum;
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
T approximateFunctionChebyshevOMP(const T& x, const std::vector<T>& coefficients) {
    int num_threads = omp_get_max_threads();
    std::vector<T> local_sums(num_threads, T(0.0));

    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        #pragma omp for
        for (size_t k = 0; k < coefficients.size(); ++k) {
            local_sums[thread_id] += coefficients[k] * chebyshevPolynomial((int)k, x);
        }
    }

    T sum = T(0.0);
    for (T local_sum : local_sums) {
        sum += local_sum;
    }

    return sum;
}

template<typename T>
std::vector<T> LoadCoeff(int maxCoefficient, std::string fileCoefficients, int numPoints, std::function<T(T)> f) {
    std::vector<T> coefficients = {};
    static std::vector<T> loadedCoefficients = {};
    if(!loadedCoefficients.empty())
        loadCoefficients(loadedCoefficients, fileCoefficients);

    if (loadedCoefficients.size() <= maxCoefficient) {
        std::cout << "Вычисление коэффициентов" << std::endl;
        coefficients.resize(maxCoefficient);
        std::copy(loadedCoefficients.begin(), loadedCoefficients.end(), coefficients.begin());
        for (int k = loadedCoefficients.size(); k < maxCoefficient; k++) {
            coefficients[k] = calculateCoefficientChebyshevOMP(k, numPoints, f);
            std::cout << "Вычислен " << k+1 << " коэффициент" << std::endl;
        }
        saveCoefficients(coefficients, fileCoefficients);
        std::cout << "Вычисление окончено" << std::endl;
    } else {
        std::cout << "Коэффициенты были загружены из файла" << std::endl;
        coefficients = loadedCoefficients;
    }
    return coefficients;
}

template<typename T>
std::vector<DataResult<T>> WorkChebyshev(T x, int maxCoefficient, int numPoints, std::function<T(T)> f, T result_x, bool isParallel, std::string fileCoefficients) {
    std::vector<T> coefficients = LoadCoeff<T>(maxCoefficient, fileCoefficients, numPoints, f);
    std::vector<DataResult<T>> results;
    std::cout << "Chebyshev:" << std::endl;
    for (int k = 0; k < maxCoefficient; k++) {
        std::vector<T> currentCoefficients(coefficients.begin(), coefficients.begin() + k + 1);

        T approxValue;
        auto start = std::chrono::high_resolution_clock::now();
        if (isParallel) {
            approxValue = approximateFunctionChebyshevOMP(x, currentCoefficients);
        } else {
            approxValue = approximateFunctionChebyshev(x, currentCoefficients);
        }
        auto end = std::chrono::high_resolution_clock::now();

        auto result = DataResult<T>(approxValue, abs(result_x - approxValue), x, k, end - start);
        results.emplace_back(result);

        std::cout << result << std::endl;
    }

    return results;
}