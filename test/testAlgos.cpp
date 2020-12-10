//
// Created by Aiday on 26/11/2020.
//

#include "../src/PowerAlgo.h"
#include "../src/InversePowerAlgo.h"
#include "../src/ShiftedPowerAlgo.h"
#include "../src/Matrix.h"
#include "../src/EigenMethods.h"
#include <utility>
#include <gtest/gtest.h>

using namespace EigenMethods;

/// Googletest fixture class to test the algorithms to compute eigenvalues
class testAlgos : public ::testing::Test {
protected:
    Matrix<double> matrix;
    int maxIter = 100;
    double tol = 0.001;
    Vector<double> initial_x0;

    virtual void SetUp() {
        // EXAMPLE 13 on following link : http://faculty.bard.edu/~belk/math213s14/Eigenvalues.pdf
        matrix.Allocate(2, 2);
        matrix(0, 0) = 2, matrix(0, 1) = 2;
        matrix(1, 0) = 1, matrix(1, 1) = 3;

        initial_x0.Allocate(2);
        initial_x0(0) = 2, initial_x0(1) = 0;
    }
};

/// Googletest class to test the power method
TEST_F(testAlgos, testPowerAlgo) {
    // EXAMPLE 13 solution
    int largest_eigenval = 4;
    Vector<double> largest_eigenvect(2);
    largest_eigenvect(0) = 1 / sqrt(2), largest_eigenvect(1) = 1 / sqrt(2);

    // Running the power algorithm
    PowerAlgo<double, double> power_algorithm_instance = PowerAlgo<double, double>(matrix, maxIter, tol);
    std::pair<double, Vector<double>> eigenval_and_vect = power_algorithm_instance.findLargest(initial_x0);
    std::cout << "Estimated Eigenvalue for the power method: " << std::get<0>(eigenval_and_vect) << std::endl;

    // Checking tolerance of eigenvalues and eigenvectors
    EXPECT_NEAR(std::get<0>(eigenval_and_vect), largest_eigenval, tol);
}

/// Googletest class to test the inverse power method
TEST_F(testAlgos, testInversePowerAlgo) {
        // EXAMPLE 13 solution
    double smallest_eigenval = 1;
    Vector<double> smallest_eigenvect(2);
    smallest_eigenvect(0) = 2/sqrt(5), smallest_eigenvect(1) = -1/sqrt(5);

   // running the inverse power algorithm
    InversePowerAlgo<double, double> inv_power_algorithm_instance = InversePowerAlgo<double, double>(matrix, maxIter, tol);
    std::pair<double, Vector<double>> eigenval_and_vect = inv_power_algorithm_instance.findSmallest(initial_x0);
    std::cout << "Estimated Eigenvalue for the inverse power method: " << std::get<0>(eigenval_and_vect) << std::endl;

    // checking tolerance of eigenvalues and eigenvectors
    EXPECT_NEAR(std::get<0>(eigenval_and_vect), smallest_eigenval, tol);
}

/// Googletest class to test the shifted inverse power method
TEST_F(testAlgos, testShiftedPowerAlgo) {
    double shift1 = 2; // eigenvalue closest will be 1
    double shift2 = 3.5; // eigenvalue closest will be 4

    // running the shifted inverse power algorithm
    ShiftedPowerAlgo<double, double> shifted_power_algorithm_instance = ShiftedPowerAlgo<double, double>(matrix, maxIter, tol);
    std::pair<double, Vector<double>> eigenval_and_vect1 = shifted_power_algorithm_instance.findEigenvalue(initial_x0, shift1);
    std::pair<double, Vector<double>> eigenval_and_vect2 = shifted_power_algorithm_instance.findEigenvalue(initial_x0, shift2);
    std::cout << "Estimated Smallest Eigenvalue for the shifted method: " << std::get<0>(eigenval_and_vect1) << std::endl;
    std::cout << "Estimated Largest Eigenvalue for the shifted method: " << std::get<0>(eigenval_and_vect2) << std::endl;

    // checking tolerance of eigenvalues and eigenvectors
    EXPECT_NEAR(std::get<0>(eigenval_and_vect1), 1, tol);
    EXPECT_NEAR(std::get<0>(eigenval_and_vect2), 4, tol);
}

/// Test divergence when eigenvalue is complex
TEST_F(testAlgos, testComplexDivergence) {
    // this matrix has complex eigenvalues,
    // with a real initial vector method can't converge
    Matrix<double> matrix2(2, 2);
    matrix2(0, 0) = 3, matrix2(0, 1) = -2;
    matrix2(1, 0) = 1, matrix2(1, 1) = 1;

    PowerAlgo<double, double> power_algorithm_instance = PowerAlgo<double, double>(matrix2, maxIter, tol);
    std::pair<double, Vector<double>> eigenval_and_vect = power_algorithm_instance.findLargest(initial_x0);

    // max iteration number should be reached
    EXPECT_EQ(maxIter, power_algorithm_instance.GetMaxIterationNumber());
}