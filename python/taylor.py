import sympy as sp
import helper

x = sp.symbols('x')

f = sp.sin(x)
x0 = sp.pi / 2
max_derivative = 100

def taylor_coefficients(func, n, xc):
    coefficients = []
    for k in range(n + 1):
        derivative = sp.diff(func, x, k)
        derivative_at_x0 = derivative.subs(x, xc)
        coefficient = derivative_at_x0 / sp.factorial(k)
        coefficients.append(coefficient)

    return coefficients

coefficients = taylor_coefficients(f, max_derivative, x0)

print("Коэффициенты ряда Тейлора для sin(x) в точке x0 = pi/2:")
helper.print_coefficients(coefficients)
