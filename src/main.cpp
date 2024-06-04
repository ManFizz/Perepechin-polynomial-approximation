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

bigfloat_t sinTaylor(bigfloat_t x, int maxTerms) {
    bigfloat_t x0 = bigfloat_t::pi / 2;
    bigfloat_t sum = 0;

    // Коэффициенты ряда Тейлора для sin(x) в точке x0 = pi/2
    bigfloat_t coefficients[] = {bigfloat_t("1"), bigfloat_t("0"), bigfloat_t("-0.5"), bigfloat_t("0"), bigfloat_t("0.0416666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666666667"), bigfloat_t("0"), bigfloat_t("-0.00138888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888888889"), bigfloat_t("0"), bigfloat_t("0.0000248015873015873015873015873015873015873015873015873015873015873015873015873015873015873015873015873015873015873015873016"), bigfloat_t("0"), bigfloat_t("-0.00000027557319223985890652557319223985890652557319223985890652557319223985890652557319223985890652557319223985890652557319224"), bigfloat_t("0"), bigfloat_t("0.00000000208767569878680989792100903212014323125434236545347656458767569878680989792100903212014323125434236545347656458767569879"), bigfloat_t("0"), bigfloat_t("-0.0000000000114707455977297247138516979786821056662326503596344866186136027405868675709945551215392485233755075024916294757564598834"), bigfloat_t("0"), bigfloat_t("0.0000000000000477947733238738529743820749111754402759693764984770275775566780857786148791439796730802021807312812603817894823185828477"), bigfloat_t("0"), bigfloat_t("-0.000000000000000156192069685862264622163643500573334235194040844696168554106791129995473461254835532941837191932291700594083275550924339"), bigfloat_t("0"), bigfloat_t("0.000000000000000000411031762331216485847799061843614037461036949591305706721333660868409140687512725086689045241927083422616008619870853523"), bigfloat_t("0"), bigfloat_t("-0.000000000000000000000889679139245057328674889744250246834331248808639189841388168097117768702786824080274218712644863816932069282726993189444"), bigfloat_t("0"), bigfloat_t("0.00000000000000000000000161173757109611834904871330480117180132472610260722797352929003101045054852685521788807737797982575531171971508513258957"), bigfloat_t("0"), bigfloat_t("-0.00000000000000000000000000247959626322479746007494354584795661742265554247265842081429235540069315157977725828934981227665500817187648474635783011"), bigfloat_t("0"), bigfloat_t("0.00000000000000000000000000000327988923706983791015204172731211192780774542655113547726758248068874755499970536810760557179451720657655619675444157422"), bigfloat_t("0"), bigfloat_t("-0.0000000000000000000000000000000037699876288159056438529215256461056641468338236219948014569913571135029367812705380547190480396749500879956284533811198")};

    for (int n = 0; n < maxTerms; ++n) {
        sum += coefficients[n] * pow(x - x0, n);
    }

    return sum;
}
void startCustomTaylor() {
    bigfloat_t x = bigfloat_t::pi / 3;
    int maxTerms = 7;
    bigfloat_t sinApprox = sinTaylor(x, maxTerms);
    std::cout << "Approximate sin = " << static_cast<long double>(sinApprox) << std::endl;
    bigfloat_t diff = abs(fsin<bigfloat_t>(x) - sinApprox);
    std::cout << "Diff = " << std::scientific << static_cast<long double>(diff) << std::endl;

}

int main() {
    bigfloat_t x = bigfloat_t("1");

    const int maxCoefficient = 14;
    const int numPoints = 50;

    start<bigfloat_t>(maxCoefficient, numPoints, x);
    startCustomTaylor();

    std::cout << "Данные сохранены в файл" << std::endl;
    return 0;
}