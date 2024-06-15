#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include "include/helper.hpp"

#include "include/Legendre.hpp"
#include "include/Chebyshev.hpp"
#include "include/Taylor.hpp"

#include "bignum.h"

void saveToFile(const char* fileName, std::vector<DataResult<bigfloat_t>>& dataResults) {
    std::ostringstream oss;
    oss << "step;result;time" << std::endl;

    for(const auto& r : dataResults) {
        oss << std::setw(2) << r.step << ";"
        << toString(r.difference, 150) << ";"
        << std::setprecision(6) << r.computationTime.count()
        << std::endl;
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

template<typename T>
T fsin(T x) {
    return sin(std::move(x), 150);
}

template<typename T>
T fcos(T x) {
    return cos(std::move(x), 150);
}

void test_cos(const int maxCoefficient, const int numPoints, bigfloat_t x, bigfloat_t result) {
    std::cout << "--  Start testing cos  --" << std::endl;
    std::cout << "Result: " << toString(result, 150) << std::endl;

    std::vector<DataResult<bigfloat_t>> dataResultsLegendre = WorkLegendre<bigfloat_t>(x, maxCoefficient, numPoints, fcos<bigfloat_t>, result, false, LegendreCosFileName);
    saveToFile("legendre_cos.csv", dataResultsLegendre);

    std::vector<DataResult<bigfloat_t>> dataResultsChebyshev = WorkChebyshev<bigfloat_t>(x, maxCoefficient, numPoints, fcos<bigfloat_t>, result, false, ChebyshevCosFileName);
    saveToFile("chebyshev_cos.csv", dataResultsChebyshev);

    std::vector<DataResult<bigfloat_t>> dataResultsTaylor = WorkTaylor(x, maxCoefficient, result, approximateCosTaylor);
    saveToFile("taylor_cos.csv", dataResultsTaylor);
}

void test_sin(const int maxCoefficient, const int numPoints, bigfloat_t x, bigfloat_t result) {
    std::cout << "--  Start testing sin  --" << std::endl;
    std::cout << "Result: " << toString(result, 150) << std::endl;

    std::vector<DataResult<bigfloat_t>> dataResultsLegendre = WorkLegendre<bigfloat_t>(x, maxCoefficient, numPoints, fsin<bigfloat_t>, result, false, LegendreSinFileName);
    saveToFile("legendre_sin.csv", dataResultsLegendre);

    std::vector<DataResult<bigfloat_t>> dataResultsChebyshev = WorkChebyshev<bigfloat_t>(x, maxCoefficient, numPoints, fsin<bigfloat_t>, result, false, ChebyshevSinFileName);
    saveToFile("chebyshev_sin.csv", dataResultsChebyshev);

    std::vector<DataResult<bigfloat_t>> dataResultsTaylor = WorkTaylor(x, maxCoefficient, result, approximateSinTaylor);
    saveToFile("taylor_sin.csv", dataResultsTaylor);
}

void test_cos_omp(const int maxCoefficient, const int numPoints, bigfloat_t x, bigfloat_t result) {
    std::cout << "--  Start testing cos omp  --" << std::endl;
    std::cout << "Result: " << toString(result, 150) << std::endl;

    std::vector<DataResult<bigfloat_t>> dataResultsLegendreOMP = WorkLegendre<bigfloat_t>(x, maxCoefficient, numPoints, fcos<bigfloat_t>, result, true, LegendreCosFileName);
    saveToFile("legendre_cos_omp.csv", dataResultsLegendreOMP);

    std::vector<DataResult<bigfloat_t>> dataResultsChebyshevOMP = WorkChebyshev<bigfloat_t>(x, maxCoefficient, numPoints, fcos<bigfloat_t>, result, true, ChebyshevCosFileName);
    saveToFile("chebyshev_cos_omp.csv", dataResultsChebyshevOMP);
}

void test_sin_omp(const int maxCoefficient, const int numPoints, bigfloat_t x, bigfloat_t result) {
    std::cout << "--  Start testing sin omp  --" << std::endl;
    std::cout << "Result: " << toString(result, 150) << std::endl;

    std::vector<DataResult<bigfloat_t>> dataResultsLegendreOMP = WorkLegendre<bigfloat_t>(x, maxCoefficient, numPoints, fsin<bigfloat_t>, result, true, LegendreSinFileName);
    saveToFile("legendre_sin_omp.csv", dataResultsLegendreOMP);

    std::vector<DataResult<bigfloat_t>> dataResultsChebyshevOMP = WorkChebyshev<bigfloat_t>(x, maxCoefficient, numPoints, fsin<bigfloat_t>, result, true, ChebyshevSinFileName);
    saveToFile("chebyshev_sin_omp.csv", dataResultsChebyshevOMP);
}

const bigfloat_t dotti = bigfloat_t("0.73908513321516064165531208767387340401341175890075746496568063577328465488354759459937610693176653184980124664398716302771490369130842031578044057462077868852490389153928943884509523480133563127677223158095635377657245120437341993643351253840978003434064670047940214347808027180188377113613820420663163350372779916967312232300613886582036217708109978970626842405880948986832618606004858989585487257367640150752276081803914595181016281591200964616460675440513264151710644662811093608258487837138395556");
const bigfloat_t result_sin_dotti = bigfloat_t("0.67361202918321483798668632649947915576381006324563461167129415074298319680434205756865987254499083621813584242002576882874918682369156027722129167141202959783080248342205558200180587832088191609565958630174311767831827545530383713935689567163198569612047443746570454427272544946788390174032078634636807560004022161336572297769140375686838943598682945839428016636173327319033980405604038938951618108883169625242256649059744788462547393420267330533052743264179917607536411059535143722416357253787025157");

int main() {
    const int maxCoefficient = 30;
    const int numPoints = 200;

    test_cos(maxCoefficient, numPoints, dotti, dotti);
    test_sin(maxCoefficient, numPoints, dotti, result_sin_dotti);

    test_cos_omp(maxCoefficient, numPoints, dotti, dotti);
    test_sin_omp(maxCoefficient, numPoints, dotti, result_sin_dotti);

    return 0;
}