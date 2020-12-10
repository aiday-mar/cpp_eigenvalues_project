//
// Created by Aiday on 26/11/2020.
//

#include "../src/EigenMethods.h"
#include "../src/Matrix.h"
#include "../src/Vector.h"
#include <gtest/gtest.h>

using namespace EigenMethods;

TEST(testEigenMethods, testRayleighQuotient) {
    Matrix<double> matrix_rayleigh(2, 2);
    matrix_rayleigh.fill(1.0);
    Vector<double> vector_rayleigh(2);
    vector_rayleigh.fill(1.0);

    EXPECT_EQ(RayleighQuotient(matrix_rayleigh, vector_rayleigh), 2);
}
