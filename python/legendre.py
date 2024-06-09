import sympy as sp
from helper import print_coefficients

def legendre_coefficients(func, n, a, b):
    x = sp.symbols('x')
    x_scaled = (2 * x - b - a) / (b - a)

    coefficients = []
    for k in range(n + 1):
        Pk = sp.legendre(k, x_scaled)
        integral = sp.integrate(func.subs(x, (b - a) / 2 * x_scaled + (b + a) / 2) * Pk * (b - a) / 2, (x, -1, 1))
        ak = (2 * k + 1) * integral
        coefficients.append(ak)

    return coefficients


x = sp.symbols('x')
f = sp.sin(x)

max_degree = 6

coefficients = legendre_coefficients(f, max_degree, 0, sp.pi / 2)

print("Коэффициенты ряда Лежандра для sin(x) на интервале [0, pi/2]:")
print_coefficients(coefficients)
