#ifndef VKR_1_DATARESULT_H
#define VKR_1_DATARESULT_H

#include <vector>
#include <chrono>

template<typename T>
class DataResult {
public:
    T result;                                           // Результат вычисления
    T difference;                                       // Отклонение
    T x;                                                // Входной параметр
    int step;                                           // Степень при которой было произведено вычисление
    std::chrono::duration<double> computationTime;      // Время вычисления

    DataResult(T res, T diff, T x, int s, std::chrono::duration<double> t)
    : result(res), difference(diff), x(x), step(s), computationTime(t) {}

    friend std::ostream& operator<<(std::ostream& os, const DataResult<T>& dr) {
        os << "Result: " << toString(dr.result, 300)
           << ", Diff: " << std::scientific << std::setprecision(0) << toString(dr.difference, 300)
           << ", Step: " << dr.step
           << ", Time: " << std::fixed << std::setprecision(6) << dr.computationTime.count() << "s";
        return os;
    }
};

#endif //VKR_1_DATARESULT_H
