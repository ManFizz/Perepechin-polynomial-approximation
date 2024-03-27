#include "../include/LegendrePolynomials.h"
#include "../include/GramSchmidtProcess.h"
#include <iostream>
#include <iomanip>

// Функция для печати полиномов
void printPolynomials(const std::vector<std::vector<double>>& polynomials) {
    for (const auto& poly : polynomials) {
        for (size_t i = 0; i < poly.size(); i++) {
            std::cout << std::fixed << std::setprecision(4) << poly[i];
            if (i != poly.size() - 1) std::cout << "x^" << i << " + ";
        }
        std::cout << "\n";
    }
}

int main() {
    const int n = 7; // Генерация и ортогонализация первых n полиномов Лежандра
    std::vector<std::vector<double>> polynomials;

    for (int i = 0; i <= n; ++i) {
        polynomials.push_back(generateLegendrePolynomial(i));
    }

    std::vector<std::vector<double>> orthogonalPolynomials = gramSchmidtProcess(polynomials);

    std::cout << "Оригинальные полиномы Лежандра:\n";
    printPolynomials(polynomials);
    std::cout << "\nОртогонализированные полиномы:\n";
    printPolynomials(orthogonalPolynomials);

    return 0;
}
