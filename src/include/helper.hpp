#ifndef VKR_HELPER_HPP
#define VKR_HELPER_HPP

#include <vector>
#include <algorithm>
#include <iostream>
#include "bignum.h"

std::string toString(const bigfloat_t& number, int precision) {
    bool isNegative = number < 0;
    bigfloat_t absNumber = isNegative ? -number : number;

    bigfloat_t scale = pow(bigfloat_t(10), precision);
    bigfloat_t scaled = absNumber * scale;
    bigfloat_t integerPart = floor(scaled);

    std::ostringstream stream;
    stream << std::fixed << std::setprecision(0) << static_cast<long double>(integerPart);
    std::string result = stream.str();

    if (precision > 0) {
        if (result.length() <= precision) {
            std::string zeros(precision + 1 - result.length(), '0');
            result = "0." + zeros + result;
        } else {
            result.insert(result.length() - precision, ".");
        }
    }

    if (isNegative) {
        result = "-" + result;
    }

    return result;
}

template<typename T>
std::vector<T> gaussianElimination(std::vector<std::vector<T>>& A, std::vector<T>& b) {
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        // Partial pivoting
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        for (int k = i; k < n; ++k) std::swap(A[maxRow][k], A[i][k]);
        std::swap(b[maxRow], b[i]);

        // Eliminate column i
        for (int k = i + 1; k < n; ++k) {
            T c = -A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                if (i == j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
            b[k] += c * b[i];
        }
    }

    // Back substitution
    std::vector<T> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i] / A[i][i];
        for (int k = i - 1; k >= 0; k--) {
            b[k] -= A[k][i] * x[i];
        }
    }

    return x;
}

// Аппроксимация методом наименьших квадратов
template<typename T, typename PolynomialFunc>
std::vector<T> fitLeastSquares(const std::vector<T>& xValues, const std::vector<T>& yValues, int degree, PolynomialFunc polyFunc) {
    int n = xValues.size();
    std::vector<std::vector<T>> A(n, std::vector<T>(degree + 1));
    std::vector<T> b(n);

    for (int i = 0; i < n; ++i) {
        b[i] = yValues[i];
        for (int j = 0; j <= degree; ++j) {
            A[i][j] = polyFunc(j, xValues[i]);
        }
    }

    std::vector<std::vector<T>> ATA(degree + 1, std::vector<T>(degree + 1, T(0)));
    std::vector<T> ATb(degree + 1, T(0));

    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k < n; ++k) {
                ATA[i][j] += A[k][i] * A[k][j];
            }
        }
        for (int k = 0; k < n; ++k) {
            ATb[i] += A[k][i] * b[k];
        }
    }

    return gaussianElimination(ATA, ATb);
}

#endif // VKR_HELPER_HPP
