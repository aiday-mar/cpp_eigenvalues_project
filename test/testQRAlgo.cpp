//
// Created by Aiday on 26/11/2020.
//

#include "../src/QRAlgo.h"
#include "../src/EigenMethods.h"
#include "../src/Matrix.h"
#include <utility>
#include <gtest/gtest.h>

using namespace EigenMethods;

class testQRAlgo : public ::testing::Test {
protected:
    // EXAMPLE on following link : http://www.math.utah.edu/~gustafso/s2019/2270/labs/lab7-polyroot-qrmethod.pdf
    Matrix<double> A;
    int max_iter = 10000;
    double tol = 0.2;
    double values_tol = 0.5;

    virtual void SetUp() {
        A.Allocate(3, 3);
        A(0,0) = 0, A(0,1) = 1, A(0,2) = 0;
        A(1,0) = 0, A(1,1) = 0, A(1,2) = 1;
        A(2,0) = -6, A(2,1) = 5, A(2,2) = 2;
    }
};

/// Googletest class to test the power method for calculating eigenvalues
TEST_F(testQRAlgo, testQRFactorization) {

    QRAlgo<double, double> qr_algorithm_instance = QRAlgo<double, double>(A, max_iter, tol);

    std::pair<Matrix<double>, Matrix<double>> qr_pair = qr_algorithm_instance.QRFactorization(A);
    std::cout << "Q" << std::endl;
    std::get<0>(qr_pair).Print();
    std::cout << std::endl;
    std::cout << "R:" << std::endl;
    std::get<1>(qr_pair).Print();
    std::cout << std::endl;

    // Check values of Q
    EXPECT_NEAR(std::get<0>(qr_pair)(0, 0), 0, values_tol);
    EXPECT_NEAR(std::get<0>(qr_pair)(0, 1), 1, values_tol);
    EXPECT_NEAR(std::get<0>(qr_pair)(0, 2), 0, values_tol);
    EXPECT_NEAR(std::get<0>(qr_pair)(1, 0), 0, values_tol);
    EXPECT_NEAR(std::get<0>(qr_pair)(1, 1), 0, values_tol);
    EXPECT_NEAR(std::get<0>(qr_pair)(1, 2), 1, values_tol);
    EXPECT_NEAR(std::get<0>(qr_pair)(2, 0), -1, values_tol);
    EXPECT_NEAR(std::get<0>(qr_pair)(2, 1), 0, values_tol);
    EXPECT_NEAR(std::get<0>(qr_pair)(2, 2), 0, values_tol);

    // Check values of R
    EXPECT_NEAR(std::get<1>(qr_pair)(0, 0), 6, values_tol);
    EXPECT_NEAR(std::get<1>(qr_pair)(0, 1), -5, values_tol);
    EXPECT_NEAR(std::get<1>(qr_pair)(0, 2), -2, values_tol);
    EXPECT_NEAR(std::get<1>(qr_pair)(1, 0), 0, values_tol);
    EXPECT_NEAR(std::get<1>(qr_pair)(1, 1),  1, values_tol);
    EXPECT_NEAR(std::get<1>(qr_pair)(1, 2),  0, values_tol);
    EXPECT_NEAR(std::get<1>(qr_pair)(2, 0),  0, values_tol);
    EXPECT_NEAR(std::get<1>(qr_pair)(2, 1),  0, values_tol);
    EXPECT_NEAR(std::get<1>(qr_pair)(2, 2),  1, values_tol);

}

TEST_F(testQRAlgo, testFindAllEigenValues) {

    QRAlgo<double, double> qr_algorithm_instance = QRAlgo<double, double>(A, max_iter, tol);
    std::pair<Vector<double>, Matrix<double>> pair = qr_algorithm_instance.ComputeEigenvalues();
    std::cout << "Vector of eigenvalues :  " << std::endl;
    std::get<0>(pair).Print();
    std::cout << "Matrix A at the end of the algorithm :  " << std::endl;
    std::get<1>(pair).Print();

    // Check eigenvalues
    EXPECT_NEAR(std::get<0>(pair)(0), 2.991, values_tol);
    EXPECT_NEAR(std::get<0>(pair)(1), -1.991, values_tol);
    EXPECT_NEAR(std::get<0>(pair)(2), 1, values_tol);
}
