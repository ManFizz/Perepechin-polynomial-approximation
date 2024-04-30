#include <cmath>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <cxxabi.h>
#include <memory>
#include <cstdlib>

#include "../include/Legendre.hpp"
#include "../include/Chebyshev.hpp"
#include "../include/Taylor.hpp"

template<typename T>
void saveToFile(const char* fileName, std::vector<DataResult<T>> dataResults) {
    std::ostringstream oss;
    oss << "step;x;result" << std::endl;

    for(const auto& result : dataResults) {
        oss << std::setw(2) << result.step << ";" << result.x << ";" << result.result << std::endl;
    }

    std::string outputString = oss.str();
    //std::cout << outputString;

    std::ofstream outputFile(fileName);
    if (outputFile.is_open()) {
        outputFile << outputString;
        outputFile.close();
        std::cout << "Данные сохранены в файл" << std::endl;
    } else {
        std::cerr << "Ошибка открытия файла" << std::endl;
    }
}

template <typename T>
std::string getTypeName() {
    int status = 0;
    std::unique_ptr<char, void(*)(void*)> res{
            abi::__cxa_demangle(typeid(T).name(), NULL, NULL, &status),
            std::free
    };
    return (status == 0) ? res.get() : typeid(T).name();
}

template<typename T>
std::string generateFileName(const std::string& methodName, const std::string& functionName) {
    std::string typeName = getTypeName<T>();
    if (typeName.find("long double") != std::string::npos) {
        typeName = "long_double";
    }
    else if (typeName.find("double") != std::string::npos) {
        typeName = "double";
    }
    return methodName + "_" + typeName + "_" + functionName + ".csv";
}

template<typename T>
T fsin(T x) {
    return std::sin(x);
}

template<typename T>
T factorial(int n) {
    T result = 1;
    for (int i = 2; i <= n; ++i) {
        result *= i;
    }
    return result;
}

template<typename T>
T sinTaylorApproximation(T x, int nTerms) {
    T sum = 0;
    for (int n = 0; n < nTerms; ++n) {
        T term = (n % 2 == 0 ? 1 : -1) * std::pow(x, 2 * n + 1) / factorial<T>(2 * n + 1);
        sum += term;
    }
    return sum;
}

template<typename T>
void start() {
    const T x = T(1);
    const int maxCoefficient = 20;
    const int numPoints = 100;
    const T result_x = fsin(x);

    //Legendre
    std::vector<DataResult<T>> dataResultsLegendre = WorkLegendre<T>(x, maxCoefficient, numPoints, fsin<T>, result_x);
    std::string fileNameLegendre = generateFileName<T>("legendre", "sin");
    saveToFile(fileNameLegendre.c_str(), dataResultsLegendre);

    //Chenyshev
    std::vector<DataResult<T>> dataResultsChebyshev = WorkChebyshev<T>(x, maxCoefficient, numPoints, fsin<T>, result_x);
    std::string fileNameChebyshev = generateFileName<T>("chebyshev", "sin");
    saveToFile(fileNameChebyshev.c_str(), dataResultsChebyshev);

    //Taylor
    auto fsinTaylor = [](T x) -> T { return sinTaylorApproximation<T>(x, 10); };
    std::vector<DataResult<T>> dataResultsTaylor = WorkTaylor<T>(x, maxCoefficient, fsinTaylor, result_x);
    std::string fileNameTaylor = generateFileName<T>("taylor", "sin");
    saveToFile(fileNameTaylor.c_str(), dataResultsTaylor);
}

int main() {
    start<long double>();
    start<double>();
    return 0;
}
