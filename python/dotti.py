import mpmath

mpmath.mp.dps = 500

dottie = mpmath.mpf(0.73908513321516064165531208767387340401341175890075746496568063577328465488354759459937610693176653184980124664398716302771490369130842031578044057462077868852490389153928943884509523480133563127677223158095635377657245120437341993643351253840978003434064670047940214347808027180188377113613820420663163350372779916967312232300613886582036217708109978970626842405880948986832618606004858989585487257367640150752276081803914595181016281591200964616460675440513264151710644662811093608258487837138395556)


def find_dottie():
    x = dottie
    for _ in range(10000):
        x = mpmath.cos(x)
    return x


dottie_number = find_dottie()
print(f"Значение косинуса числа дотти с высокой точностью: {dottie_number}")

sin_dottie = mpmath.sin(dottie_number)
print(f"Значение sin в точке числа дотти: {sin_dottie}")
