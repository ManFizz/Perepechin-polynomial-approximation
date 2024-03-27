#include "../include/GramSchmidtProcess.h"
#include "../include/LegendrePolynomials.h"
#include <vector>
#include <cmath>

// Функция для скалярного произведения векторов
double dotProduct(const std::vector<double>& v1, const std::vector<double>& v2) {
    double sum = 0.0;
    for (size_t i = 0; i < v1.size() && i < v2.size(); i++) {
        sum += v1[i] * v2[i];
    }
    return sum;
}

// Функция для умножения вектора на скаляр
std::vector<double> scalarMultiply(const std::vector<double>& v, double scalar) {
    std::vector<double> result(v.size(), 0.0);
    for (size_t i = 0; i < v.size(); i++) {
        result[i] = v[i] * scalar;
    }
    return result;
}

// Функция для вычитания двух векторов
std::vector<double> vectorSubtract(const std::vector<double>& v1, const std::vector<double>& v2) {
    std::vector<double> result(v1.size(), 0.0);
    for (size_t i = 0; i < v1.size() && i < v2.size(); i++) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

//Функция процесса Грама-Шмидта
std::vector<std::vector<double>> gramSchmidtProcess(const std::vector<std::vector<double>>& polynomials) {
    std::vector<std::vector<double>> orthogonalPolynomials;
    for (const auto& p : polynomials) {
        std::vector<double> u = p;
        for (const auto& q : orthogonalPolynomials) {
            double dotQQ = dotProduct(q, q);
            if (std::fabs(dotQQ) > std::numeric_limits<double>::epsilon()) {
                u = vectorSubtract(u, scalarMultiply(q, dotProduct(p, q) / dotQQ));
            }
        }
        orthogonalPolynomials.push_back(u);
    }
    return orthogonalPolynomials;
}