
#include "include/helper.hpp"

#include "include/Legendre.hpp"
#include "include/Chebyshev.hpp"
#include "include/Taylor.hpp"

#include "params.hpp"

void test_legendre_cos(bigfloat_t& x, bigfloat_t& r) {
    auto dr = WorkLegendre<bigfloat_t>(x, fcos<bigfloat_t>, r, false, LegendreCosFileName);
    saveToFile("legendre_cos.csv", dr);
}

void test_chebyshev_cos(bigfloat_t& x, bigfloat_t& r) {
    auto dr = WorkChebyshev<bigfloat_t>(x, fcos<bigfloat_t>, r, false, ChebyshevCosFileName);
    saveToFile("chebyshev_cos.csv", dr);
}

void test_taylor_cos(bigfloat_t& x, bigfloat_t& r) {
    auto dr = WorkTaylor(x, maxCoefficient, r, approximateCosTaylor);
    saveToFile("taylor_cos.csv", dr);
}

void test_legendre_sin(bigfloat_t& x, bigfloat_t& r) {
    auto dr = WorkLegendre<bigfloat_t>(x, fsin<bigfloat_t>, r, false, LegendreSinFileName);
    saveToFile("legendre_sin.csv", dr);
}

void test_chebyshev_sin(bigfloat_t& x, bigfloat_t& r) {
    auto dr = WorkChebyshev<bigfloat_t>(x, fsin<bigfloat_t>, r, false, ChebyshevSinFileName);
    saveToFile("chebyshev_sin.csv", dr);
}

void test_taylor_sin(bigfloat_t x, bigfloat_t r) {
    auto dr = WorkTaylor(x, maxCoefficient, r, approximateSinTaylor);
    saveToFile("taylor_sin.csv", dr);
}



void test_legendre_cos_omp(bigfloat_t& x, bigfloat_t& r) {
    auto dr = WorkLegendre<bigfloat_t>(x, fcos<bigfloat_t>, r, true, LegendreCosFileName);
    saveToFile("legendre_cos_omp.csv", dr);
}

void test_chebyshev_cos_omp(bigfloat_t& x, bigfloat_t& r) {
    auto dr = WorkChebyshev<bigfloat_t>(x, fcos<bigfloat_t>, r, true, ChebyshevCosFileName);
    saveToFile("chebyshev_cos_omp.csv", dr);
}

void test_legendre_sin_omp(bigfloat_t& x, bigfloat_t& r) {
    auto dr = WorkLegendre<bigfloat_t>(x, fsin<bigfloat_t>, r, true, LegendreSinFileName);
    saveToFile("legendre_sin_omp.csv", dr);
}

void test_chebyshev_sin_omp(bigfloat_t& x, bigfloat_t& r) {
    auto dr = WorkChebyshev<bigfloat_t>(x, fsin<bigfloat_t>, r, true, ChebyshevSinFileName);
    saveToFile("chebyshev_sin_omp.csv", dr);
}

void test_cos(bigfloat_t& x, bigfloat_t& r) {
    std::cout << "--  Start testing cos  --" << std::endl;
    test_legendre_cos(x, r);
    test_chebyshev_cos(x, r);
    test_taylor_cos(x, r);
    std::cout << "--  End testing cos  --" << std::endl;
}

void test_sin(bigfloat_t& x, bigfloat_t& r) {
    std::cout << "--  Start testing sin  --" << std::endl;
    test_legendre_sin(x, r);
    test_chebyshev_sin(x, r);
    test_taylor_sin(x, r);
    std::cout << "--  End testing sin  --" << std::endl;
}

void test_cos_omp(bigfloat_t& x, bigfloat_t& r) {
    std::cout << "--  Start testing cos omp  --" << std::endl;
    test_legendre_cos_omp(x, r);
    test_chebyshev_cos_omp(x, r);
    std::cout << "--  End testing cos omp  --" << std::endl;
}


void test_sin_omp(bigfloat_t& x, bigfloat_t& r) {
    std::cout << "--  Start testing sin omp  --" << std::endl;
    test_legendre_sin_omp(x, r);
    test_chebyshev_sin_omp(x, r);
    std::cout << "--  End testing sin omp  --" << std::endl;
}

void fastTest() {

    std::cout << "sin omp legendre: " << toString(sin(dotti, 250), 250) << std::endl;

    {
        std::vector<bigfloat_t> coefficients = {};
        loadCoefficients(coefficients, "python_legendre_coefficients_sin.txt");
        std::vector<bigfloat_t> currentCoefficients(coefficients.begin(), coefficients.begin() + 60);
        auto value = approximateFunctionLegendreOMP(dotti, currentCoefficients);
        std::cout << "sin omp legendre: " << toString(abs(value)) << std::endl;
    }
}

void compareCoefficients() {
    const int checkCoefficient = 60;
    {
        std::vector<bigfloat_t> coefficients = {};
        loadCoefficients(coefficients, "legendre_coefficients_cos.txt");
        std::vector<bigfloat_t> currentCoefficients(coefficients.begin(), coefficients.begin() + checkCoefficient + 1);
        auto value = approximateFunctionLegendreOMP(dotti, currentCoefficients);
        std::cout << "c++ coeff cos_omp legendre: " << countZeros(abs(value - dotti)) << std::endl;
    }

    {
        std::vector<bigfloat_t> coefficients = {};
        loadCoefficients(coefficients, "python_legendre_coefficients_cos.txt");
        std::vector<bigfloat_t> currentCoefficients(coefficients.begin(), coefficients.begin() + checkCoefficient + 1);
        auto value = approximateFunctionLegendreOMP(dotti, currentCoefficients);
        std::cout << "python coeff cos_omp legendre: " << countZeros(abs(value - dotti)) << std::endl;
    }

    {
        std::vector<bigfloat_t> coefficients = {};
        loadCoefficients(coefficients, "chebyshev_coefficients_cos.txt");
        std::vector<bigfloat_t> currentCoefficients(coefficients.begin(), coefficients.begin() + checkCoefficient + 1);
        auto value = approximateFunctionChebyshevOMP(dotti, currentCoefficients);
        std::cout << "c++ coeff cos_omp cheb: " << countZeros(abs(value - dotti)) << std::endl;
    }

    {
        std::vector<bigfloat_t> coefficients = {};
        loadCoefficients(coefficients, "python_chebyshev_coefficients_cos.txt");
        std::vector<bigfloat_t> currentCoefficients(coefficients.begin(), coefficients.begin() + checkCoefficient + 1);
        auto value = approximateFunctionChebyshevOMP(dotti, currentCoefficients);
        std::cout << "python coeff cos_omp cheb: " << countZeros(abs(value - dotti)) << std::endl;
    }
}