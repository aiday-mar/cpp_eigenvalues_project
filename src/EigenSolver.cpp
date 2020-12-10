//
// Created by Aiday on 23/11/2020.
//
#include "EigenSolver.h"

template<typename T, typename S>
EigenSolver<T,S>::EigenSolver(const Matrix<S>& matrix, int maxIter, T tol) : mMatrix(matrix), mMaxIter(maxIter), mTol(tol) {
}

// only implementing get methods as the default constructor
// does create mMatrix, mTol and mMaxIter
template<typename T, typename S>
Matrix<S> EigenSolver<T,S>::GetMatrix() {
    return mMatrix;
}

template<typename T, typename S>
double EigenSolver<T,S>::GetTolerance() const {
    return mTol;
}

template<typename T, typename S>
int EigenSolver<T,S>::GetMaxIterationNumber() const{
    return mMaxIter;
}

// supported types (double, float)
template class EigenSolver<double, double>;
template class EigenSolver<float, float>;