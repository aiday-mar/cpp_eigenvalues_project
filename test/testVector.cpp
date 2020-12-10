//
// Created by Carmine on 29/11/2020.
//

#include "../src/Vector.h"
#include "../src/Exception.h"
#include <gtest/gtest.h>
#include <complex>

/// Googletest fixture class to test the vector class
class testVector : public ::testing::Test {
protected:
    Vector<double> vector1;
    Vector<double> vector2;
    Vector<double> vector3;
    Vector<std::complex<double>> vector1_c;
    Vector<std::complex<double>> vector2_c;
    std::complex<double> a;

    virtual void SetUp() {
        vector1.Allocate(2);
        vector2.Allocate(2);
        vector3.Allocate(2);
        vector1_c.Allocate(3);
        vector2_c.Allocate(3);
        a = std::complex<double>(3. + 4.j);
    }
};

TEST_F(testVector, testGetSetEntries) {
    EXPECT_EQ(vector1(0), 0);
    EXPECT_EQ(vector1(1), 0);

    vector1.fill(5.5);
    EXPECT_EQ(vector1(0), 5.5);
    EXPECT_EQ(vector1(1), 5.5);
}

// == operator test is needed for EXPECT_EQ function for the other tests
TEST_F(testVector, testOperatorComparison) {
    vector1.fill(1.0);
    vector2.fill(1.0);
    EXPECT_EQ(vector1, vector2);
}

TEST_F(testVector, testInvalidConstructor) {
    EXPECT_THROW(Vector<double> s(-2), NegativeNumberException);
    EXPECT_THROW(Vector<double> t(-2), NegativeNumberException);
}

TEST_F(testVector, testCopyConstructor) {
    Vector<double> v(vector1);
    EXPECT_EQ(v, vector1);
}

TEST_F(testVector, testGetSize) {
    EXPECT_EQ(vector1.GetSize(), 2);
}

TEST_F(testVector, testOperatorEqual) {
    vector1.fill(1.0);
    vector2 = vector1;
    EXPECT_EQ(vector2, vector1);
}

TEST_F(testVector, testOperatorPlus) {
    vector1.fill(1.0);
    vector2.fill(1.0);
    vector3.fill(2.0);
    EXPECT_EQ(vector1 + vector2, vector3);
}

TEST_F(testVector, tematrixstOperatorMinus) {
    vector1.fill(1.0);
    vector2.fill(1.0);
    vector3.fill(0.0);
    EXPECT_EQ(vector2 - vector1, vector3);
}

TEST_F(testVector, testOperatorMinusUnary) {
    vector1.fill(-2.0);
    vector2.fill(2.0);
    EXPECT_EQ(vector1, - vector2);
}

TEST_F(testVector, testOperatorScalarMulti) {
    vector1.fill(1.0);
    vector2.fill(2.0);
    EXPECT_EQ(vector2, vector1 * 2);
}

TEST_F(testVector, testNorm) {
    vector1(0) = 8, vector1(1) = 6;
    EXPECT_EQ(vector1.Norm(), 10.0);
}

TEST_F(testVector, testNormalized) {
    vector1(0) = 8, vector1(1) = 0;
    vector2(0) = 1, vector2(1) = 0;
    vector1.Normalize();
    EXPECT_EQ(vector1, vector2);
}

TEST_F(testVector, testNormComplex) {
    std::complex<double> res = sqrt(75) + 0j;
    vector1_c.fill(a);
    EXPECT_EQ(vector1_c.Norm(), res);
}

TEST_F(testVector, testVectorComplex) {
    vector1_c.fill(a);
    Vector<std::complex<double>> res(vector1_c);
    EXPECT_EQ(vector1_c, res);
}

TEST_F(testVector, testNormalizeAndScalarMultiplicationComplex) {
    vector1_c.fill(a);
    Vector<std::complex<double>> res(vector1_c);
    std::complex<double> var = 1./vector1_c.Norm();
    res = res * (var);
    vector1_c.Normalize();
    EXPECT_EQ(vector1_c, res);
}

TEST_F(testVector, testUnaryMinusComplex) {
    vector1_c.fill(a);
    vector2_c = - vector1_c;
    vector1_c.fill(-a);
    EXPECT_EQ(vector1_c, vector2_c);
}

TEST_F(testVector, testPlusComplex) {
    vector1_c.fill(a);
    vector2_c.fill(2.*a);
    vector1_c = vector2_c + vector1_c;
    vector2_c.fill(3.*a);
    EXPECT_EQ(vector1_c, vector2_c);
}

TEST_F(testVector, testMinusComplex) {
    vector1_c.fill(a);
    vector2_c.fill(2.*a);
    vector1_c = vector2_c - vector1_c;
    vector2_c.fill(a);
    EXPECT_EQ(vector1_c, vector2_c);
}

TEST_F(testVector, testScalarProductComplex) {
    vector1_c.fill(a);
    std::complex<double> b;
    b = vector1_c * vector1_c;
    std::complex<double> res(-21.+72.j);
    EXPECT_EQ(b, res);
}


