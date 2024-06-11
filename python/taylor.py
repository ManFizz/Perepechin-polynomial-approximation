import sympy as sp
from helper import print_coefficients

def taylor_coefficients_sin(func, n, xc):
    x = sp.symbols('x')
    coefficients = []
    for k in range(n + 1):
        derivative = sp.diff(func, x, k)
        derivative_at_x0 = derivative.subs(x, xc)
        coefficient = derivative_at_x0 / sp.factorial(k)
        coefficients.append(coefficient)

    return coefficients

x = sp.symbols('x')
f = sp.sin(x)
x0 = 0.5
max_derivative = 100

coefficients = taylor_coefficients_sin(f, max_derivative, x0)

print("Коэффициенты ряда Тейлора для sin(x) в точке x0 = pi/2:")
print_coefficients(coefficients)
