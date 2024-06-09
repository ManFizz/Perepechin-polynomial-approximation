import sympy as sp
from helper import print_coefficients

def chebyshev_coefficients(func, n, a=0, b=sp.pi / 2):
    x = sp.symbols('x')
    x_j = [sp.cos(sp.pi * (j + 0.5) / n) for j in range(n)]
    x_j_transformed = [0.5 * (xj + 1) * (b - a) + a for xj in x_j]

    f_j = [func.subs(x, xjt) for xjt in x_j_transformed]

    coefficients = []
    for k in range(n):
        term = (2 / n) * sum(fj * sp.cos(sp.pi * k * (j + 0.5) / n) for j, fj in enumerate(f_j))
        coefficients.append(term)

    coefficients[0] /= 2

    return coefficients

x = sp.symbols('x')
f = sp.sin(x)

max_derivative = 10

coefficients = chebyshev_coefficients(f, max_derivative, 0, sp.pi / 2)

print("Коэффициенты ряда Чебышева для sin(x) в точке x0 = pi/2:")
print_coefficients(coefficients)
