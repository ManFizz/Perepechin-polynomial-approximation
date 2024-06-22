from mpmath import mp, cos, sin
from decimal import Decimal, getcontext

mp.dps = 500


def calculate_coefficient_legendre(k, num_points, f):
    def legendre_polynomial(x):
        return mp.legenp(k, 0, x)

    def integrand(x):
        return f(x) * legendre_polynomial(x)

    coefficient = mp.quad(integrand, [-1, 1], method='gauss-legendre', maxdegree=num_points)

    coefficient *= (2 * k + 1) / 2

    return coefficient


def near_zero(value, tolerance=1e-150):
    return abs(value) < tolerance


def save_coefficients_to_file(coefficients, filename):
    getcontext().prec = 150
    with open(filename, 'w') as file:
        for coef in coefficients:
            if near_zero(coef):
                file.write("0.0\n")
            else:
                file.write(format(Decimal(str(coef)), 'f') + '\n')


def evaluate_legendre_series(coefficients, x):
    mp.dps = 150
    result = mp.mpf(0)

    for k, c_k in enumerate(coefficients):
        P_k = mp.legenp(k, 0, x)
        result += c_k * P_k

    return result

num_points = 130
max_degree = 100

sin_coefficients = []
for k in range(max_degree + 1):
    coef = calculate_coefficient_legendre(k, num_points, sin)
    sin_coefficients.append(coef)

x = 0.73908513321516064165531208767387340401341175890075746496568063577328465488354759459937610693176653184980124664398716302771490369130842031578044057462077868852490389153928943884509523480133563127677223158095635377657245120437341993643351253840978003434064670047940214347808027180188377113613820420663163350372779916967312232300613886582036217708109978970626842405880948986832618606004858989585487257367640150752276081803914595181016281591200964616460675440513264151710644662811093608258487837138395556
# save_coefficients_to_file(sin_coefficients, '../data/python_legendre_coefficients_sin.txt')

value = evaluate_legendre_series(sin_coefficients, x)
print(format(Decimal(str(value)), 'f'))


# cos_coefficients = []
# for k in range(max_degree + 1):
#     coef = calculate_coefficient_legendre(k, num_points, cos)
#     cos_coefficients.append(coef)
# save_coefficients_to_file(cos_coefficients, '../data/python_legendre_coefficients_cos.txt')
#
# print("Коэффициенты сохранены в файлы ../data/python_legendre_coefficients_sin.txt и ../data/python_legendre_coefficients_cos.txt")
