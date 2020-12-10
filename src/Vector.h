//
// Created by carmine on 29/11/2020.
//

#ifndef PSC_PROJECT_EIGENVALUES_VECTOR_H
#define PSC_PROJECT_EIGENVALUES_VECTOR_H

#include <iostream>
#include <cmath>

template<typename T> class Vector {
public:
    // Public constants and variables
    static double TOL;
    // Constructors
    Vector();
    Vector(const int size);
    Vector(const Vector<T>& otherVector);
    // Destructor
    virtual ~Vector();
    // Getters
    int GetSize() const;
    // Setters
    // nothing
    // Methods
    void Allocate(const int size);
    void Print(std::ostream &output = std::cout) const;
    T* begin();
    T* end();
    T Norm() const;
    void Normalize();
    void fill(T value);
    // Operator
    T& operator()(const int i);
    Vector<T>& operator=(const Vector<T> &v);
    bool operator==(const Vector<T> &v) const;
    bool operator!=(const Vector<T> &v) const;
    Vector<T> operator-() const;
    Vector<T> operator+(const Vector<T> &v) const;
    Vector<T> operator-(const Vector<T> &v) const;
    Vector<T> operator*(const T a) const;
    T operator*(const Vector<T> &v) const;
//    Vector<T> operator*(Matrix<T> &M) const;
private:
    T* mData;
    int mSize;
    bool mAllocated;
};

#endif //PSC_PROJECT_EIGENVALUES_VECTOR_H
