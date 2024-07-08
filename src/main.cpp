#include "params.hpp"
#include "tests.cpp"

int main() {
    fastTest();
    compareCoefficients();
    test_legendre_cos_omp(dotti, dotti);
    TestApproxChebCos();

    return 0;
}