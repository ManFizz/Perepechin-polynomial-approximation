#include "../include/LegendrePolynomials.h"
#include <vector>
#include <cmath>

// Функция для вычисления полиномов Лежандра
std::vector<double> generateLegendrePolynomial(int n) {
    std::vector<double> coefficients(n + 1, 0.0);

    if (n == 0) {
        coefficients[0] = 1;
        return coefficients;
    }

    if (n == 1) {
        coefficients[1] = 1;
        return coefficients;
    }

    std::vector<double> Pn_1 = generateLegendrePolynomial(n - 1);
    std::vector<double> Pn_2 = generateLegendrePolynomial(n - 2);

    // Рекурсивное соотношение для вычисления полиномов Лежандра
    for (int i = 0; i <= n - 1; ++i) {
        coefficients[i + 1] += ((2.0 * n - 1) / n) * Pn_1[i];
        coefficients[i] -= ((n - 1.0) / n) * Pn_2[i];
    }

    return coefficients;
}
