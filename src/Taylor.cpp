#include <vector>

#include "include/DataResult.h"
#include "include/helper.hpp"


bigfloat_t approximateSinTaylor(bigfloat_t& x, int terms) {
    static bigfloat_t sum = 0;
    static int last_n = 0;
    static bigfloat_t last_term = 0;

    if (terms == 0) {
        sum = 0;
        last_n = 0;
        last_term = x;
        sum += last_term;
    } else {
        for (int n = last_n + 1; n <= terms; ++n) {
            last_term *= -1 * x * x / ((2 * n) * (2 * n + 1));
            sum += last_term;
        }
    }
    last_n = terms;
    return sum;
}

bigfloat_t approximateCosTaylor(bigfloat_t& x, int terms) {
    static bigfloat_t sum = 0;
    static int last_n = 0;
    static bigfloat_t last_term = 0;

    if (terms == 0) {
        sum = 0;
        last_n = 0;
        last_term = 1;
        sum += last_term;
    } else {
        for (int n = last_n + 1; n <= terms; ++n) {
            last_term *= -1 * x * x / ((2 * n - 1) * (2 * n));
            sum += last_term;
        }
    }
    last_n = terms;
    return sum;
}

std::vector<DataResult<bigfloat_t>> WorkTaylor(bigfloat_t& x, int maxCoefficient, bigfloat_t& result_x, std::function<bigfloat_t(bigfloat_t&, int)> f) {
    std::vector<DataResult<bigfloat_t>> results;

    std::cout << "Work Taylor" << std::endl;
    for (int k = 0; k < maxCoefficient; ++k) {

        auto start = std::chrono::high_resolution_clock::now();
        bigfloat_t approxValue = f(x, k);
        auto end = std::chrono::high_resolution_clock::now();

        auto result = DataResult<bigfloat_t>(approxValue, abs(result_x - approxValue), x, k, end - start);
        results.emplace_back(result);

        PRINT_IF_LOGS_ENABLED(result)
    }
    f(x, -1);
    PRINT_END()

    return results;
}