#ifndef VKR_1_DATARESULT_H
#define VKR_1_DATARESULT_H

#include <vector>

template<typename T>
class DataResult {
public:
    T result;
    T x;
    int step;

    DataResult(T res, T val, int s) : result(res), x(val), step(s) {}

    static void AddData(std::vector<DataResult<T>>& data, T result, T x, int step) {
        data.emplace_back(result, x, step);
    }
};

#endif //VKR_1_DATARESULT_H
