//
// Created by carmine on 24/11/2020.
//

#include "Matrix.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include "Vector.h"
#include "Exception.h"
#include <complex>

// Constants
// To check if two numbers are the same
template<typename T>
double Matrix<T>::TOL = 1e-8;

// Private methods
// Algorithm that computes the determinant of a matrix (recursive)
template<typename T>
T CalculateDeterminant(Matrix<T>& A, const int A_size) {
    T det = 0.;

    if (A_size == 2) {
        det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
        return det;
    } else {
        for (int i = 0; i < A_size; i++) {
            Matrix<T> A_reduced(A_size-1, A_size-1);
            A.ReduceMatrix(A_reduced, 0, i);
            det += A(0, i) * std::pow(-1., i) * CalculateDeterminant(A_reduced, A_size-1);
        }
        return det;
    };
}

/// Constructors/Destructors
// Constructor for initialization
template<typename T>
Matrix<T>::Matrix() { mAllocated = false; }

// To allocate memory and initialise entries to zero
template<typename T>
Matrix<T>::Matrix(const int numRows, const int numCols) {
    mAllocated = false;
    Allocate(numRows, numCols);
}

// Copy constructor
// To allocate for a new memory and copy the entries into this matrix from the other one
template<typename T>
Matrix<T>::Matrix(const Matrix<T> &otherMatrix) {
    mNumRows = otherMatrix.mNumRows;
    mNumCols = otherMatrix.mNumCols;

    mData = new T[mNumRows * mNumCols];
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            mData[i * mNumCols + j] = otherMatrix.mData[i * mNumCols + j];
        }
    }
    mAllocated = true;
}

// Load matrix from a text file (complex type not supported)
template<typename T>
Matrix<T>::Matrix(const std::string fileName) {
    // Check if file exists, the file has to be placed in the src folder
    std::ifstream read_file(fileName.c_str());
    if (read_file.is_open() == false) {
        throw FileNotOpenException();
    }

    // Read first line with number of rows and columns
    int numRows;
    int numCols;
    read_file >> numRows >> numCols;

    // Allocate when matrix is empty
    if (!mAllocated) {
        Allocate(numRows, numCols);
    }

    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumRows; j++) {
            read_file >> mData[i * mNumCols + j];
        }
    }

    read_file.close();

    std::cout << fileName << " read successfully" << std::endl;
}

// Destructor
// To free memory
template<typename T>
Matrix<T>::~Matrix() {
    delete[] mData;
}

/// Getters
// Get the number of rows of the matrix
template<typename T>
int Matrix<T>::GetNumberOfRows() const {
    return mNumRows;
}

// Get the number of cols of the matrix
template<typename T>
int Matrix<T>::GetNumberOfCols() const {
    return mNumCols;
}

// Get the column j of the matrix
template <typename T>
Vector<T> Matrix<T>::getColumn(const int j) const {
    Vector<T> result(GetNumberOfRows());
    for (int i = 0; i < GetNumberOfRows(); i++) {
        result(i) = mData[i*GetNumberOfCols()+j];
    }
    return result;
}

template <typename T>
Vector<T> Matrix<T>::getRow(const int i) const {
    Vector<T> result(GetNumberOfCols());
    for (int j = 0; j < GetNumberOfCols(); j++) {
        result(j) = mData[i*GetNumberOfCols()+j];
    }
    return result;
}

template<typename T>
T* Matrix<T>::GetData() const {
    return mData;
}

/// Setters
// nothing

/// Methods
// Allocate memory and fill the matrix if the instance of the class was created using the constructor Matrix()
template<typename T>
void Matrix<T>::Allocate(const int numRows, const int numCols) {
    if(mAllocated) {
        throw AlreadyAllocatedException();
    }
    if (numRows < 1 || numCols < 1) {
        throw NegativeNumberException();
    }

    mNumRows = numRows;
    mNumCols = numCols;
    mData = new T [mNumRows*mNumCols];
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            mData[i*mNumCols+j] = 0;
        }
    }

    mAllocated = true;
}

// Print the matrix
template<typename T>
void Matrix<T>::Print(std::ostream &output) const {
    for (int i = 0; i < GetNumberOfRows(); i++) {
        for (int j = 0; j < GetNumberOfCols(); j++) {
            output << mData[i*GetNumberOfCols()+j] << " ";
        }
        output << std::endl;
    }
}

// Return the pointer of the first entry of the matrix
template<typename T>
T *Matrix<T>::begin() {
    return &mData[0];
}

// Return the pointer of the last entry of the matrix
template<typename T>
T *Matrix<T>::end() {
    return &mData[0] + GetNumberOfCols()*GetNumberOfRows();
}

// Compute the reduced matrix
template<typename T>
void Matrix<T>::ReduceMatrix(Matrix<T>& A_reduced, const int row_skipped, const int col_skipped) const {
    if (GetNumberOfRows() < 1 || GetNumberOfCols() < 1) {
        throw WrongDimensionException();
    }
    if (GetNumberOfRows() != (A_reduced.GetNumberOfRows()+1) || GetNumberOfCols() != (A_reduced.GetNumberOfCols()+1)) {
        throw WrongDimensionException();
    }

    for (int row = 0; row < row_skipped; row++) {
        for (int col = 0; col < col_skipped; col++) {
            A_reduced(row, col) = mData[row*GetNumberOfCols()+col];
        }
        for(int col = col_skipped+1; col < GetNumberOfCols(); col++) {
            A_reduced(row, col-1) = mData[row*GetNumberOfCols()+col];
        }
    }
    for (int row = row_skipped+1; row < GetNumberOfRows(); row++) {
        for (int col = 0; col < col_skipped; col++) {
            A_reduced(row-1, col) = mData[row*GetNumberOfCols()+col];
        }
        for(int col = col_skipped+1; col < GetNumberOfCols(); col++) {
            A_reduced(row-1, col-1) = mData[row*GetNumberOfCols()+col];
        }
    }
}

// Return the transpose
template<typename T>
Matrix<T> Matrix<T>::Transpose() const {
    Matrix<T> res(mNumCols, mNumRows);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            res(j, i) = mData[i*mNumCols+j];
        }
    }

    return res;
}

// Return the determinant
template<typename T>
T Matrix<T>::Determinant() const {
    if(GetNumberOfCols() != GetNumberOfRows()) {
        throw NonSquareMatrixException();
    }

    Matrix<T> A(*this);
    T det = CalculateDeterminant(A, A.GetNumberOfRows());

    return det;
}

// Return the inverse
template<typename T>
Matrix<T> Matrix<T>::Inverse() const {
    if(GetNumberOfCols() != GetNumberOfRows()) {
        throw NonSquareMatrixException();
    }
    if(IsDiagonal()) {
        throw SingularMatrixException();
    }

    Matrix<T> res(mNumRows, mNumCols);
    const T det = Determinant();

    if (mNumRows == 2) {
        res(0, 0) = mData[3]/det;
        res(0, 1) = -mData[1]/det;
        res(1, 0) = -mData[2]/det;
        res(1, 1) = mData[0]/det;

    } else {
        Matrix<T> tmp_cofactor(mNumRows, mNumCols);
        tmp_cofactor = Cofactor();

        res = tmp_cofactor.Transpose() * (1./det);
    }

    return res;
}

// Return true if the matrix is singular
template<typename T>
bool Matrix<T>::IsSingular() const {
    if (std::abs(this->Determinant()) < TOL) {
        return true;
    } else {
        return false;
    }
}

// Return true if the matrix is symmetric
template<typename T>
bool Matrix<T>::IsSymmetric() const {
    if(GetNumberOfCols() != GetNumberOfRows()) {
        throw NonSquareMatrixException();
    }

    for (int i = 0; i < mNumRows-1; i++) {
        for (int j = i+1; j < mNumCols; j++) {
            if (std::abs(mData[i*mNumCols+j]-mData[j*mNumCols+i])  > TOL) { return false; };
        }
    }

    return true;
}

// Return true if the matrix is diagonal
template<typename T>
bool Matrix<T>::IsDiagonal() const {
    if(GetNumberOfCols() != GetNumberOfRows()) {
        throw NonSquareMatrixException();
    }

    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            if ((i != j) && (std::abs(mData[i*mNumCols+j]) > TOL)) { return false; };
        }
    }

    return true;
}

// Compute the cofactor matrix
template<typename T>
Matrix<T> Matrix<T>::Cofactor() const {
    if(GetNumberOfCols() != GetNumberOfRows()) {
        throw NonSquareMatrixException();
    }

    Matrix<T> res(mNumRows, mNumCols);

    // Looping for each element of the matrix
    for (int row = 0; row < mNumRows; row++) {
        for (int col = 0; col < mNumCols; col++) {
            Matrix<T> tmp(mNumRows-1, mNumCols-1);
            ReduceMatrix(tmp, row, col);
            res(row, col) = tmp.Determinant() * pow(-1, row+col);
        }
    }
    return res;
}

// Function which fills all entries of a vector with the same value
template<typename T>
void Matrix<T>::fill(T value) {
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            mData[i*mNumCols+j] = value;
        }
    }
}

/// Parser functions
// Load matrix from command line
template <typename T>
void Matrix<T>::InputFromCommandLine() {
    // Command line dialog to input matrix
    std::cout << "********************************** " << std::endl
              << "Input values for Matrix(" << mNumRows << "," << mNumCols<< "):" << std::endl;
    // Nested for loop to store values for new matrix
    for (int i = 0; i < GetNumberOfRows(); i++) {
        for (int j = 0; j < GetNumberOfCols(); j++) {
            std::cout << "Input [" << i << ", " << j << "]: ";
            std::cin >> mData[i*mNumCols+j];
        }
    }
    std::cout << "**********************************" << std::endl;
}

/// Operators
// Return the entry (i, j) (possible to set or get the entry)
// Overloading the round brackets
template<typename T>
T& Matrix<T>::operator()(const int i, const int j) {
    if(i < 0 || j < 0 || i >= mNumRows || j >= mNumCols) {
        throw std::out_of_range("Index out of range!");
    }

    return mData[i*mNumCols+j];
}

// Overloading the assign operator
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& M) {
    if(mNumRows != M.GetNumberOfRows() || mNumCols != M.GetNumberOfCols()) {
        throw WrongDimensionException();
    }

    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            mData[i*mNumCols+j] = M.mData[i*mNumCols+j];
        }
    }

    return *this;
}

// Overloading the == operator
template<typename T>
bool Matrix<T>::operator==(const Matrix<T> &M) const {
    if(mNumRows != M.GetNumberOfRows() || mNumCols != M.GetNumberOfCols()) {
        throw WrongDimensionException();
    }

    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            if (std::abs(mData[i*mNumCols+j]-M.mData[i*mNumCols+j]) > TOL) {
                return false;
            }
        }
    }
    return true;
}

// Overloading the != operator
template<typename T>
bool Matrix<T>::operator!=(const Matrix<T> &M) const {
    if(mNumRows != M.GetNumberOfRows() || mNumCols != M.GetNumberOfCols()) {
        throw WrongDimensionException();
    }

    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            if (std::abs(mData[i*mNumCols+j]-M.mData[i*mNumCols+j]) > TOL) {
                return true;
            }
        }
    }
    return false;
}

// Overloading the unary - operator
template<typename T>
Matrix<T> Matrix<T>::operator-() const {
    Matrix<T> res(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            res(i, j) = -mData[i*mNumCols+j];
        }
    }

    return res;
}

// Overloading the binary + operator
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &M) const {
    if(mNumRows != M.GetNumberOfRows() || mNumCols != M.GetNumberOfCols()) {
        throw WrongDimensionException();
    }

    Matrix<T> res(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            res(i, j) = mData[i*mNumCols+j] + M.mData[i*mNumCols+j];
        }
    }

    return res;
}

// Overloading the binary - operator
template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &M) const {
    if(mNumRows != M.GetNumberOfRows() || mNumCols != M.GetNumberOfCols()) {
        throw WrongDimensionException();
    }

    Matrix<T> res(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            res(i, j) = mData[i*mNumCols+j] - M.mData[i*mNumCols+j];
        }
    }

    return res;
}

// Overloading the binary * operator (scalar)
template<typename T>
Matrix<T> Matrix<T>::operator*(const T a) const {
    Matrix<T> res(mNumRows, mNumCols);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            res(i, j) = mData[i*mNumCols+j] * a;
        }
    }

    return res;
}

// Overloading  the binary * operator (matrix * vector)
template<typename T>
Vector<T> Matrix<T>::operator* (Vector<T>& v) const {
    if(v.GetSize() != mNumCols) {
        throw WrongDimensionException();
    }

    Vector<T> res(mNumRows);
    for (int i = 0; i < mNumRows; i++) {
        for (int j = 0; j < mNumCols; j++) {
            res(i) += mData[i*mNumCols+j] * v(j);
        }
    }

    return res;
}

// Overloading the binary * operator (matrix * matrix)
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& M) const{
    if(mNumCols != M.GetNumberOfRows()) {
        throw WrongDimensionException();
    }

    int numberColumnsOtherMatrix = M.GetNumberOfCols();
    // matrix of size the outer dimension
    Matrix<T> res(mNumRows, numberColumnsOtherMatrix);
    for(int i=0; i < mNumRows; i++) {
        for(int j=0; j<numberColumnsOtherMatrix; j++){
            res(i,j) = 0;
        }
    }

    for(int i = 0; i < mNumRows; i++) {
        for(int j = 0; j < numberColumnsOtherMatrix; j++){
            for(int k = 0; k < mNumCols; k++){
                res(i, j) += mData[i*mNumCols+k]*M.mData[k*mNumCols+j];
            }
        }
    }
    return res;
}

// Overloading  the binary * operator (vector * matrix)
template <typename T> Vector<T> operator*(Vector<T> &v, Matrix<T> &M) {
    if (v.GetSize() != M.GetNumberOfRows()) {
        throw WrongDimensionException();
    }

    Vector<T> res(M.GetNumberOfCols());
    for (int i = 0; i < M.GetNumberOfCols(); i++) {
        for (int j = 0; j < M.GetNumberOfRows(); j++) {
            res(i) +=  v(j) * M(j, i);
        }
    }

    return res;
}

// Prototypes (supporting double, float)
//template Vector<int> operator*(Vector<int> &v, Matrix<int> &M);
template Vector<float> operator*(Vector<float> &v, Matrix<float> &M);
template Vector<double> operator*(Vector<double> &v, Matrix<double> &M);
template Vector<std::complex<double>> operator*(Vector<std::complex<double>> &v, Matrix<std::complex<double>> &M);

// Prototypes (supporting double, float)
//template class Matrix<int>;
template class Matrix<double>;
template class Matrix<float>;
template class Matrix<std::complex<double>>;