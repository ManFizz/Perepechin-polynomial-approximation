import sympy as sp

x = sp.symbols('x')

f = sp.sin(x)
x0 = sp.pi / 2
max_derivative = 30
precision = 120

coefficients = []
for k in range(max_derivative + 1):
    derivative = sp.diff(f, x, k)
    derivative_at_x0 = derivative.subs(x, x0)
    coefficient = derivative_at_x0 / sp.factorial(k)
    coefficients.append(coefficient)

formatted_coefficients = []
for coeff in coefficients:
    decimal_coeff = str(coeff.evalf(precision))
    decimal_coeff = decimal_coeff.rstrip('0').rstrip('.') if '.' in decimal_coeff else decimal_coeff
    formatted_coefficients.append(f'bigfloat_t("{decimal_coeff}")')

print("Коэффициенты ряда Тейлора для sin(x) в точке x0 = pi/2:")
cpp_output = "bigfloat_t coefficients[] = {" + ", ".join(formatted_coefficients) + "};"
print(cpp_output)
