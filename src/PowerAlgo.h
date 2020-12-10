//
// Created by Aiday on 26/11/2020.
//

#ifndef PSC_PROJECT_EIGENVALUES_POWERALGO_H
#define PSC_PROJECT_EIGENVALUES_POWERALGO_H

#include "EigenSolver.h"
#include "EigenMethods.h"
#include "Matrix.h"
#include "Vector.h"

using namespace EigenMethods;

template<typename T, typename S> class PowerAlgo : public EigenSolver<T,S> {
public :
    /// default constructor
    PowerAlgo(const Matrix<S>& matrix, int maxIter, T tol);
    /**
     * @param x_initial
     * @return the largest eigenvalue of the matrix and the corresponding eigenvector where you use the initial estimate x_initial
     */
    std::pair<T, Vector<T>> findLargest(const Vector<T>& x_initial);
};

#endif //PSC_PROJECT_EIGENVALUES_POWERALGO_H