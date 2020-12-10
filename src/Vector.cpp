//
// Created by carmine on 29/11/2020.
//

#include "Vector.h"
#include <iostream>
#include <cmath>
#include "Exception.h"
#include <complex>

// Constants
// To check if two numbers are the same
template<typename T>
double Vector<T>::TOL = 1e-8;

/// Constructors/Destructors
// Constructor for initialization
template<typename T>
Vector<T>::Vector() { mAllocated = false; }

// To allocate memory and initialise entries to zero
template<typename T>
Vector<T>::Vector(const int size) {
    mAllocated = false;
    Allocate(size);
}

// Copy constructor
// To allocate for a new memory and copy the entries into this matrix from the other one
template<typename T>
Vector<T>::Vector(const Vector<T>& otherVector) {
    mSize = otherVector.mSize;
    mData = new T [mSize];
    for (int i = 0; i < mSize; i++) {
        mData[i] = otherVector.mData[i];
    }

    mAllocated = true;
}

// Destructor
// To free memory
template<typename T>
Vector<T>::~Vector() {
    delete[] mData;
}

/// Getters
// Get the number of rows of the matrix
template<typename T>
int Vector<T>::GetSize() const {
    return mSize;
}

/// Setters
// nothing

/// Methods
// Allocate memory and fill the vector if the instance of the class was created using the constructor Vector()
template<typename T>
void Vector<T>::Allocate(const int size) {
    if (mAllocated) {
        throw AlreadyAllocatedException();
    }
    if (size < 1) {
        throw NegativeNumberException();
    }

    if (!mAllocated) {
        mSize = size;
        mData = new T [mSize];
        for (int i = 0; i < mSize; i++) {
            mData[i] = 0;
        }

        mAllocated = true;
    }
}

// Print the matrix
template<typename T>
void Vector<T>::Print(std::ostream &output) const {
    for (int i = 0; i < mSize; i++) {
        output << mData[i] << std::endl;
    }
}

// Return the pointer of the first entry of the matrix
template<typename T>
T* Vector<T>::begin() {
    return &mData[0];
}


// Return the pointer of the last entry of the matrix
template<typename T>
T* Vector<T>::end() {
    return &mData[0] + mSize;
}

// Return the norm
template<typename T>
T Vector<T>::Norm() const {
    T res = 0;
    for (int i = 0; i < mSize; i++) {
        res += pow(std::abs(mData[i]), 2);
    }
    return sqrt(res);
}

// Change the vector turning it into a normalized one
template <typename T>
void Vector<T>::Normalize() {
    T norm = Norm();
    for (int i = 0; i < mSize; i++) {
        mData[i] = mData[i]/norm;
    }
}

// Function which fills all entries of a Vector with the same value
template <typename T>
void Vector<T>::fill(T value) {
for (int i = 0; i < GetSize(); i++) {
    mData[i] = value;
    }
}

template<typename T>
T &Vector<T>::operator()(const int i) {
    if(i < 0 || i >= mSize) {
        throw std::out_of_range("Index out of range!");
    }

    return mData[i];
}

template<typename T>
Vector<T> &Vector<T>::operator=(const Vector<T> &v) {
    if(mSize != v.GetSize()) {
        throw WrongDimensionException();
    }

    for (int i = 0; i < mSize; i++) {
        mData[i] = v.mData[i];
    }

    return *this;
}

template<typename T>
bool Vector<T>::operator==(const Vector<T> &v) const {
    if(mSize != v.GetSize()) {
        throw WrongDimensionException();
    }

    for (int i = 0; i < mSize; i++) {
        if (std::abs(mData[i]-v.mData[i]) > TOL) {
            return false;
        }
    }
    return true;
}

template<typename T>
bool Vector<T>::operator!=(const Vector<T> &v) const {
    if(mSize != v.GetSize()) {
        throw WrongDimensionException();
    }

    for (int i = 0; i < mSize; i++) {
        if (std::abs(mData[i]-v.mData[i]) > TOL) {
            return true;
        }
    }
    return false;
}

template<typename T>
Vector<T> Vector<T>::operator-() const {
    Vector<T> res(mSize);
    for (int i = 0; i < mSize; i++) {
        res(i) = -mData[i];
    }

    return res;
}

template<typename T>
Vector<T> Vector<T>::operator+(const Vector<T> &v) const {
    if(mSize != v.GetSize()) {
        throw WrongDimensionException();
    }

    Vector<T> res(mSize);
    for (int i = 0; i < mSize; i++) {
        res(i) = mData[i] + v.mData[i];
    }

    return res;
}

template<typename T>
Vector<T> Vector<T>::operator-(const Vector<T> &v) const {
    if(mSize != v.GetSize()) {
        throw WrongDimensionException();
    }

    Vector<T> res(mSize);
    for (int i = 0; i < mSize; i++) {
        res(i) = mData[i] - v.mData[i];
    }

    return res;
}

template<typename T>
Vector<T> Vector<T>::operator*(const T a) const {
    Vector<T> res(mSize);
    for (int i = 0; i < mSize; i++) {
        res(i) = mData[i] * a;
    }

    return res;
}

// Scalar product
template<typename T>
T Vector<T>::operator*(const Vector<T> &v) const {
    T res = 0.;
    for (int i = 0; i < mSize; i++) {
        res += mData[i] * v.mData[i];
    }

    return res;
}

// Prototypes (supporting int, double, float)
//template class Vector<int>;
template class Vector<double>;
template class Vector<float>;
template class Vector<std::complex<double>>;


