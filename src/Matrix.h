//
// Created by carmine on 24/11/2020.
//

#ifndef INC_8_PROJECT_MATRIX_H
#define INC_8_PROJECT_MATRIX_H

#include <iostream>
#include <cmath>
#include <vector>
#include "Vector.h"

// Matrix of general type T (int, float, double)
template<typename T> class Matrix {
public:
    // Public constants and variables
    static double TOL;
    // Constructors
    Matrix();
    Matrix(const int numRows,const int numCols);
    Matrix(const Matrix<T>& otherMatrix);
    Matrix(const std::string fileName);
    // Destructor
    virtual ~Matrix();
    // Getters
    int GetNumberOfRows() const;
    int GetNumberOfCols() const;
    T* GetData() const;
    Vector<T> getColumn(const int j) const;
    Vector<T> getRow(const int i) const;
    // Setters
    // nothing
    // Methods
    void Allocate(const int numRows, const int numCols);
    void Print(std::ostream &output = std::cout) const;
    T* begin();
    T* end();
    void ReduceMatrix(Matrix<T>& A_reduced, const int row_skipped, const int col_skipped) const;
    Matrix<T> Inverse() const;
    Matrix<T> Transpose() const;
    T Determinant() const;
    bool IsSingular() const;
    bool IsSymmetric() const;
    bool IsDiagonal() const;
    Matrix<T> Cofactor() const;
    void fill(T value);
    void InputFromCommandLine();

    // Operator
    T& operator()(const int i, const int j);
    Matrix<T>& operator=(const Matrix<T>& M);
    bool operator==(const Matrix<T>& M) const;
    bool operator!=(const Matrix<T>& M) const;
    Matrix<T> operator-() const;
    Matrix<T> operator+(const Matrix<T>& M) const;
    Matrix<T> operator-(const Matrix<T>& M) const;
    Matrix<T> operator*(const T a) const;
    Vector<T> operator*(Vector<T>& v) const;
    Matrix<T> operator*(const Matrix<T>& M) const;
private:
    // Private constants and variables
    int mNumRows, mNumCols; // dimensions
    bool mAllocated;
    T* mData; // entries of the matrix
};

// declaration of the operation vector * matrix
template <typename T>
Vector<T> operator*(Vector<T> &v, Matrix<T> &M);

#endif //INC_8_PROJECT_MATRIX_H