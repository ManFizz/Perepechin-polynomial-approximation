#include <vector>
#include <cmath>
#include <functional>

#include "DataResult.h"

template<typename T>
std::vector<DataResult<T>> WorkTaylor(const T x, const int maxCoefficient, std::function<T(T)> f, const T result_x) {
    std::vector<DataResult<T>> results;
    for(int k = 1; k <= maxCoefficient; k++) {
        T result = f(x);
        DataResult<T>::AddData(results, std::abs(result_x - result), x, k);
    }
    return results;
}