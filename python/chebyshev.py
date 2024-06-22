from mpmath import mp, cos, sin
from decimal import Decimal, getcontext

mp.dps = 200


def calculate_coefficient_chebyshev(k, num_points, f):

    def chebyshev_polynomial(x):
        return mp.chebyt(k, x)

    nodes = [mp.cos(mp.pi * (i + 0.5) / num_points) for i in range(num_points)]
    weights = [mp.pi / num_points] * num_points

    coefficient = 0
    for x, w in zip(nodes, weights):
        coefficient += w * f(x) * chebyshev_polynomial(x)

    coefficient *= 2 / mp.pi

    if k == 0:
        coefficient /= 2

    return coefficient


def near_zero(value, tolerance=1e-200):
    return abs(value) < tolerance


def save_coefficients_to_file(coefficients, filename):
    getcontext().prec = 200
    with open(filename, 'w') as file:
        for coef in coefficients:
            if near_zero(coef):
                file.write("0.0\n")
            else:
                file.write(format(Decimal(str(coef)), 'f') + '\n')

num_points = 150
max_degree = 100

sin_coefficients = []
for k in range(max_degree + 1):
    coef = calculate_coefficient_chebyshev(k, num_points, sin)
    sin_coefficients.append(coef)
save_coefficients_to_file(sin_coefficients, '../data/python_chebyshev_coefficients_sin.txt')

cos_coefficients = []
for k in range(max_degree + 1):
    coef = calculate_coefficient_chebyshev(k, num_points, cos)
    cos_coefficients.append(coef)
save_coefficients_to_file(cos_coefficients, '../data/python_chebyshev_coefficients_cos.txt')

print("Коэффициенты сохранены в файлы ../data/python_chebyshev_coefficients_sin.txt и ../data/python_chebyshev_coefficients_cos.txt")