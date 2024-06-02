#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <cxxabi.h>
#include <utility>

#include "include/Legendre.hpp"
#include "include/Chebyshev.hpp"
#include "include/Taylor.hpp"

#include "bignum.h"

void saveToFile(const char* fileName, std::vector<DataResult<bigfloat_t>>& dataResults) {
    std::ostringstream oss;
    oss << "step;result" << std::endl;

    for(const auto& result : dataResults) {
        oss << std::setw(2) << result.step << ";" << static_cast<long double>(result.result) << std::endl;
    }

    std::string outputString = oss.str();

    std::ofstream outputFile(fileName);
    if (outputFile.is_open()) {
        outputFile << outputString;
        std::cout << outputString; //For debug
        outputFile.close();
    } else {
        std::cerr << "Ошибка открытия файла" << std::endl;
    }
}

bigfloat_t complexFunction(bigfloat_t x) {
    if (x <= 0) {
        return sin(5 * x, 54);
    } else {
        return 1 / (1 + 25 * x * x);
    }
}

template<typename T>
T fsin(T x) {
    //return complexFunction(x);
    return sin(std::move(x), 54);
}

template<typename T>
void start(const int maxCoefficient, const int numPoints, T x) {
    T result_x = fsin(x);

    //Legendre
    std::vector<DataResult<T>> dataResultsLegendre = WorkLegendre<T>(x, maxCoefficient, numPoints, fsin<T>, result_x);
    saveToFile("legendre_sin.csv", dataResultsLegendre);

    //Chebyshev
    std::vector<DataResult<T>> dataResultsChebyshev = WorkChebyshev<T>(x, maxCoefficient, numPoints, fsin<T>, result_x);
    saveToFile("chebyshev_sin.csv", dataResultsChebyshev);

    //Taylor
    std::vector<DataResult<T>> dataResultsTaylor = WorkTaylor<T>(x, maxCoefficient, fsin<T>, result_x);
    saveToFile("taylor_sin.csv", dataResultsTaylor);
}

int main() {
    bigfloat_t x = bigfloat_t("0.765");

    const int maxCoefficient = 21;
    const int numPoints = 100;

    start<bigfloat_t>(maxCoefficient, numPoints, x);

    std::cout << "Данные сохранены в файл" << std::endl;
    return 0;
}