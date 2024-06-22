#include "params.hpp"
#include "tests.cpp"

int main() {
    //fastTest();
    //compareCoefficients();
    //test_legendre_cos_omp(dotti, dotti);

    test_cos(dotti, dotti);
    test_sin(dotti, result_sin_dotti);

    test_cos_omp(dotti, dotti);
    test_sin_omp(dotti, result_sin_dotti);

    return 0;
}