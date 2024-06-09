#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>

#include "include/Legendre.hpp"
#include "include/Chebyshev.hpp"
#include "include/Taylor.hpp"

#include "bignum.h"

void saveToFile(const char* fileName, std::vector<DataResult<bigfloat_t>>& dataResults) {
    std::ostringstream oss;
    oss << "step;result" << std::endl;

    for(const auto& result : dataResults) {
        oss << std::setw(2) << result.step << ";" << static_cast<long double>(result.result) << std::endl;
        std::cout << result << std::endl; //For debug
    }

    std::string outputString = oss.str();

    std::ofstream outputFile(fileName);
    if (outputFile.is_open()) {
        outputFile << outputString;
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

void start(const int maxCoefficient, const int numPoints, bigfloat_t x) {
    bigfloat_t result_x = fsin(x);

    //Legendre
    //std::vector<DataResult<T>> dataResultsLegendre = WorkLegendre<bigfloat_t>(x, maxCoefficient, numPoints, fsin<T>, result_x);
    //saveToFile("legendre_sin.csv", dataResultsLegendre);

    //Chebyshev
    //std::vector<DataResult<T>> dataResultsChebyshev = WorkChebyshev<bigfloat_t>(x, maxCoefficient, numPoints, fsin<T>, result_x);
    //saveToFile("chebyshev_sin.csv", dataResultsChebyshev);

    //Taylor
    std::vector<DataResult<bigfloat_t>> dataResultsTaylor = WorkTaylor(x, maxCoefficient, fsin<bigfloat_t>, result_x);
    saveToFile("taylor_sin.csv", dataResultsTaylor);
}

const bigfloat_t dotti = bigfloat_t("0.73908513321516064165531208767387340401341175890075746496568063577328465488354759459937610693176653184980124664398716302771490369130842031578044057462077868852490389153928950778242826380949334414333694");

int main() {
    bigfloat_t x = dotti;

    const int maxCoefficient = 14;
    const int numPoints = 50;
    start(maxCoefficient, numPoints, x);

    return 0;
}