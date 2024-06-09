precision = 200

def print_coefficients(coefficients):
    formatted_coefficients = []
    for coeff in coefficients:
        decimal_coeff = coeff.evalf(precision)
        formatted_str = f"{decimal_coeff:.{precision}f}"
        formatted_str = formatted_str.rstrip('0').rstrip('.')
        formatted_coefficients.append(f'bigfloat_t("{formatted_str}")')

    cpp_output = "bigfloat_t coefficients[] = {" + ", ".join(formatted_coefficients) + "};"
    print(cpp_output)

