#ifndef VKR_HELPER_HPP
#define VKR_HELPER_HPP

#include <vector>
#include <algorithm>
#include <iostream>
#include "bignum.h"
#include <omp.h>
#include "DataResult.h"

#include <fstream>
#include <sstream>

std::string toString(const bigfloat_t& number, int precision = TO_STRING_PRECISION) {
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
            std::string zeros(precision - result.length(), '0');
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

int countZeros(const bigfloat_t& val) {
    auto number = toString(val);
    size_t pointIndex = number.find('.');
    if (pointIndex == std::string::npos) {
        return 0;
    }

    int zeroCount = 0;
    for (size_t i = pointIndex + 1; i < number.length(); ++i) {
        if (number[i] == '0') {
            zeroCount++;
        } else {
            break;
        }
    }

    return zeroCount;
}

template<typename T>
std::vector<T> gaussianElimination(std::vector<std::vector<T>>& A, std::vector<T>& b) {
    int n = A.size();

    for (int i = 0; i < n; ++i) {
        // Partial pivoting
        int maxRow = i;
        #pragma omp parallel for shared(maxRow)
        for (int k = i + 1; k < n; ++k) {
            #pragma omp critical
            {
                if (abs(A[k][i]) > abs(A[maxRow][i])) {
                    maxRow = k;
                }
            }
        }

        for (int k = i; k < n; ++k) std::swap(A[maxRow][k], A[i][k]);
        std::swap(b[maxRow], b[i]);

        // Eliminate column i
        #pragma omp parallel for
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
        #pragma omp parallel for
        for (int k = 0; k < i; ++k) {
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

    #pragma omp parallel for collapse(2)
    for (int i = 0; i <= degree; ++i) {
        for (int j = 0; j <= degree; ++j) {
            for (int k = 0; k < n; ++k) {
                #pragma omp critical
                ATA[i][j] += A[k][i] * A[k][j];
            }
        }
    }

    std::vector<T> ATb(degree + 1, T(0));
    #pragma omp parallel for
    for (int i = 0; i <= degree; ++i) {
        for (int k = 0; k < n; ++k) {
            #pragma omp critical
            ATb[i] += A[k][i] * b[k];
        }
    }

    return gaussianElimination(ATA, ATb);
}

void createFileIfNotExists(const std::string& workFile) {
    std::ifstream infile(workFile);
    if (!infile.good()) {
        std::ofstream file(workFile);
        if (!file.is_open())
            std::cerr << "Не удалось создать файл: " << workFile << std::endl;
        file.close();
    }
}

template<typename T>
void saveCoefficients(const std::vector<T>& coefficients, const std::string& filename) {
    const std::string workPath = PathToData + filename;
    createFileIfNotExists(workPath);

    std::ofstream outFile(workPath);
    if (!outFile) {
        std::cerr << "Не удалось открыть файл с коэфициентами." << std::endl;
        return;
    }

    for (const auto& coef : coefficients) {
        outFile << toString(coef, 250) << std::endl;
    }
    outFile.close();
}

template<typename T>
bool loadCoefficients(std::vector<T>& coefficients, const std::string& filename) {
    const std::string workPath = PathToData + filename;
    createFileIfNotExists(workPath);

    std::ifstream inFile(workPath);
    if (!inFile) {
        std::cerr << "Не удалось открыть файл с коэфициентами." << std::endl;
        return false;
    }

    std::string value;
    while (inFile >> value) {
        coefficients.push_back(T(value));
    }
    inFile.close();
    return true;
}

void saveToFile(const char* fileName, std::vector<DataResult<bigfloat_t>>& dataResults) {
const std::string workPath = PathToData + fileName;

std::ostringstream oss;
oss << "step;result;difference;time;x" << std::endl;

for(const auto& r : dataResults) {
oss << std::setw(2) << r.step << ";"
<< toString(r.result) << ";"
<< toString(r.difference) << ";"
<< std::setprecision(6) << r.computationTime.count() << ";"
<< toString(r.x)
<< std::endl;
}

std::string outputString = oss.str();

std::ofstream outputFile(workPath);
if (outputFile.is_open()) {
outputFile << outputString;
outputFile.close();
} else {
std::cerr << "Ошибка открытия файла" << std::endl;
}
}

#endif // VKR_HELPER_HPP
