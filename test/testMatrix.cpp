//
// Created by Carmine on 24/11/2020.
//

#include "../src/Matrix.h"
#include "../src/Vector.h"
#include "../src/Exception.h"
#include <complex>
#include <fstream>
#include <gtest/gtest.h>

/// Googletest fixture class to test the matrix class
class testMatrix : public ::testing::Test {
protected:
    Matrix<double> matrix1;
    Matrix<double> matrix2;
    Matrix<double> matrix3;
    Matrix<std::complex<double>> matrix1_c;
    Matrix<std::complex<double>> matrix2_c;
    Matrix<std::complex<double>> matrix3_c;
    Matrix<std::complex<double>> matrix4_c;
    Matrix<std::complex<double>> matrix5_c;
    std::complex<double> a;

    Vector<double> vector1;
    Vector<double> vector2;
    Vector<std::complex<double>> vector1_c;
    Vector<std::complex<double>> vector2_c;

    virtual void SetUp() {
        matrix1.Allocate(2, 2);
        matrix2.Allocate(2, 2);
        matrix3.Allocate(2, 2);
        matrix1_c.Allocate(2, 2);
        matrix2_c.Allocate(2, 2);
        matrix3_c.Allocate(2, 3);
        matrix4_c.Allocate(2, 3);
        matrix5_c.Allocate(3, 3);
        a = std::complex<double>(3. + 4.j);

        vector1.Allocate(2);
        vector2.Allocate(2);
        vector1_c.Allocate(2);
        vector2_c.Allocate(3);
    }
};

TEST_F(testMatrix, testGetSetEntries) {
    EXPECT_EQ(matrix1(0, 0), 0);
    EXPECT_EQ(matrix1(0, 1), 0);
    EXPECT_EQ(matrix1(1, 0), 0);
    EXPECT_EQ(matrix1(1, 1), 0);

    matrix1.fill(5.5);
    EXPECT_EQ(matrix1(0, 0), 5.5);
    EXPECT_EQ(matrix1(0, 1), 5.5);
    EXPECT_EQ(matrix1(1, 0), 5.5);
    EXPECT_EQ(matrix1(1, 1), 5.5);
}

// == operator test is need for EXPECT_EQ function for the other tests
TEST_F(testMatrix, testOperatorComparison) {
    matrix1.fill(1.0);
    matrix2.fill(1.0);
    EXPECT_EQ(matrix1, matrix2);
}

TEST_F(testMatrix, testInvalidConstructor) {
    EXPECT_THROW(Matrix<double> A(2, -2), NegativeNumberException);
    EXPECT_THROW(Matrix<double> B(-2, -2), NegativeNumberException);
}

TEST_F(testMatrix, testCopyConstructor) {
    Matrix<double> A(matrix1);
    EXPECT_EQ(A, matrix1);
}

TEST_F(testMatrix, testGetNumRowsAndCols) {
    EXPECT_EQ(matrix1.GetNumberOfCols(), 2);
    EXPECT_EQ(matrix1.GetNumberOfRows(), 2);
}

TEST_F(testMatrix, testOperatorEqual) {
    matrix1.fill(1.0);
    matrix2 = matrix1;
    EXPECT_EQ(matrix2, matrix1);
}

TEST_F(testMatrix, testOperatorPlus) {
    matrix1.fill(1.0);
    matrix2.fill(1.0);
    matrix3.fill(2.0);
    EXPECT_EQ(matrix1 + matrix2, matrix3);
}

TEST_F(testMatrix, testOperatorMinus) {
    matrix1.fill(1.0);
    matrix2.fill(1.0);
    matrix3.fill(0.0);
    EXPECT_EQ(matrix2 - matrix1, matrix3);
}

TEST_F(testMatrix, testOperatorMinusUnary) {
    matrix1.fill(-2.0);
    matrix2.fill(2.0);
    EXPECT_EQ(matrix1, - matrix2);
}

TEST_F(testMatrix, testOperatorScalarMulti) {
    matrix1.fill(2.0);
    matrix2.fill(1.0);
    EXPECT_EQ(matrix1, matrix2 * 2);
}

TEST_F(testMatrix, testOperatorMultiplication) {
    matrix1.fill(1.0);
    matrix2.fill(1.0);
    matrix3.fill(2.0);
    EXPECT_EQ(matrix1 * matrix2, matrix3);
}

TEST_F(testMatrix, testOperatorMatrixVectorMulti) {
    vector1.fill(1.0);
    matrix1.fill(1.0);
    vector2.fill(2.0);
    EXPECT_EQ(vector2, matrix1 * vector1);
}

TEST_F(testMatrix, testDeterminant) {
    matrix1.fill(1.0);
    matrix1(0, 0) = 3;
    EXPECT_EQ(matrix1.Determinant(), 2);
}

TEST_F(testMatrix, testInverse) {
    matrix1(0, 0) = 5;
    matrix1(0, 1) = 2;
    matrix1(1, 0) = -7;
    matrix1(1, 1) = -3;
    matrix2(0, 0) = 3;
    matrix2(0, 1) = 2;
    matrix2(1, 0) = -7;
    matrix2(1, 1) = -5;
    EXPECT_EQ(matrix1.Inverse(), matrix2);
}

TEST_F(testMatrix, testTranspose) {
    matrix1.fill(1.0);
    matrix1(0, 1) = 2;
    matrix2.fill(1.0);
    matrix2(1, 0) = 2;
    EXPECT_EQ(matrix1.Transpose(), matrix2);
}

TEST_F(testMatrix, testSingular) {
    matrix1.fill(1.0);
    matrix2.fill(1.0);
    matrix2(1, 0) = 0;
    EXPECT_TRUE(matrix1.IsSingular());
    EXPECT_FALSE(matrix2.IsSingular());
}

TEST_F(testMatrix, testDiagonal) {
    matrix1.fill(1.0);
    matrix1(1, 0) = 0;
    matrix1(0, 1) = 0;
    matrix2.fill(1.0);
    EXPECT_TRUE(matrix1.IsDiagonal());
    EXPECT_FALSE(matrix2.IsDiagonal());
}

TEST_F(testMatrix, testVectorComplex) {
    matrix1_c.fill(a);
    Matrix<std::complex<double>> res(matrix1_c);
    EXPECT_EQ(matrix1_c, res);
}

TEST_F(testMatrix, testSingularComplex) {
    matrix1_c.fill(a);
    matrix1_c(0, 0) = std::complex<double>(0 + 0j);
    matrix2_c.fill(a);
    EXPECT_TRUE(matrix2_c.IsSingular());
    EXPECT_FALSE(matrix1_c.IsSingular());
}

TEST_F(testMatrix, testUnaryMinusComplex) {
    matrix1_c.fill(a);
    matrix2_c = - matrix1_c;
    matrix1_c.fill(-a);
    EXPECT_EQ(matrix1_c, matrix2_c);
}

TEST_F(testMatrix, testPlusComplex) {
    matrix1_c.fill(a);
    matrix2_c.fill(2.*a);
    matrix1_c = matrix2_c + matrix1_c;
    matrix2_c.fill(3.*a);
    EXPECT_EQ(matrix1_c, matrix2_c);
}

TEST_F(testMatrix, testMinusComplex) {
    matrix1_c.fill(a);
    matrix2_c.fill(2.*a);
    matrix1_c = matrix2_c - matrix1_c;
    matrix2_c.fill(a);
    EXPECT_EQ(matrix1_c, matrix2_c);
}

TEST_F(testMatrix, testScalarMultiplicationComplex) {
    matrix1_c.fill(2.*a);
    matrix2_c.fill(a);
    matrix2_c= matrix2_c * 2.;
    EXPECT_EQ(matrix1_c, matrix2_c);
}

TEST_F(testMatrix, testMatrixTimesVectorComplex) {
    vector2_c.fill(a);
    matrix3_c.fill(a);
    vector1_c.fill(std::complex<double>(-21. + 72.j));
    EXPECT_EQ(matrix3_c * vector2_c, vector1_c);
}

TEST_F(testMatrix, testMatrixTimesMatrixComplex) {
    matrix1_c.fill(a);
    matrix3_c.fill(std::complex<double>(1. + 0j));
    matrix4_c.fill(2.*a);
    EXPECT_EQ(matrix1_c * matrix3_c, matrix4_c);
}

TEST_F(testMatrix, testVectorTimesMatrixComplex) {
    vector1_c.fill(a);
    matrix3_c.fill(a);
    vector2_c.fill(std::complex<double>(-14. + 48.j));
    EXPECT_EQ(vector1_c * matrix3_c, vector2_c);
}

TEST_F(testMatrix, testReducedMatrixComplex) {
    matrix5_c.fill(a);
    matrix5_c(0, 2) = std::complex<double>(0. + 0.j);
    matrix5_c(1, 2) = std::complex<double>(0. + 0.j);
    matrix5_c(2, 0) = std::complex<double>(0. + 0.j);
    matrix5_c(2, 1) = std::complex<double>(0. + 0.j);
    matrix5_c(1, 0) = std::complex<double>(1. + 0.j);
    matrix2_c.fill(std::complex<double>(0 + 0j));
    matrix2_c(0, 0) = std::complex<double>(1. + 0j);
    matrix2_c(1, 1) = a;
    matrix5_c.ReduceMatrix(matrix1_c, 0, 1);
    EXPECT_EQ(matrix1_c, matrix2_c);
}

TEST_F(testMatrix, testDeterminantComplex) {
    matrix5_c.fill(a);
    matrix5_c(0, 2) = std::complex<double>(0. + 0.j);
    matrix5_c(1, 2) = std::complex<double>(0. + 0.j);
    matrix5_c(2, 0) = std::complex<double>(0. + 0.j);
    matrix5_c(2, 1) = std::complex<double>(0. + 0.j);
    matrix5_c(1, 0) = std::complex<double>(1. + 0.j);
    std::complex<double> det;
    det = matrix5_c.Determinant();
    EXPECT_EQ(det, std::complex<double>(-110. + 20.j));
}

TEST_F(testMatrix, testInverseComplex) {
    matrix1_c(0, 0) = a;
    matrix1_c(0, 1) = std::complex<double>(0 + 0j);
    matrix1_c(1, 0) = std::complex<double>(1. + 0j);
    matrix1_c(1, 1) = std::complex<double>(3. - 4.j);
    matrix2_c(0, 0) = std::complex<double>(3. - 4.j);
    matrix2_c(0, 1) = std::complex<double>(0 + 0j);
    matrix2_c(1, 0) = std::complex<double>(-1. + 0j);
    matrix2_c(1, 1) = a;
    matrix2_c = matrix2_c * (1./std::complex<double>(25.+0j));
    EXPECT_EQ(matrix1_c.Inverse(), matrix2_c);
}

/// Googletest for constructor from valid text file
TEST(testMatrixLoaders, testValidFile) {
    Matrix<double> A(3, 3);
    A.fill(1.0);

    // Store Matrix in file (assuming this works)
    std::ofstream new_file;
    new_file.open("matrix.txt");

    new_file << A.GetNumberOfRows() << " " << A.GetNumberOfCols() << std::endl;

    for (int i = 0; i < A.GetNumberOfRows(); i++) {
        for (int j = 0; j < A.GetNumberOfCols(); j++) {
            new_file << A(i, j) << " ";
        }
        new_file << std::endl;
    }
    new_file.close();

    // Matrix Reader instance to store values
    Matrix<double> B("matrix.txt");

    EXPECT_EQ(A, B);

    remove( "matrix.txt");
}

/// Googletest for constructor from invalid text file
TEST(testMatrixLoaders, testInvalidFile) {
    Matrix<double> C(3, 3);
    C.fill(1.0);

    // Open empty file
    std::ofstream new_file;
    new_file.open("invalid_matrix.txt");

    // If file is not correctly will allocate 0 for missing values
    try {
        Matrix<double> D("invalid_matrix.txt");
        remove("invalid_matrix.txt");
        FAIL();
    }
    // Not specifying rows and columns should throw exception
    catch(NegativeNumberException) {
        remove("invalid_matrix.txt");
    }
}

/// Googletest for constructor from not existing text file
TEST(testMatrixLoaders, testFileNotOpenException) {
    try {
        Matrix<double> E("not_a_filename.txt");
        FAIL();
    }

    catch(FileNotOpenException) {}
}
