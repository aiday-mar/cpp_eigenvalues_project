//
// Created by Aiday on 23/11/2020.
//

#ifndef PSC_PROJECT_EIGENVALUES_EIGENSOLVER_H
#define PSC_PROJECT_EIGENVALUES_EIGENSOLVER_H

#include "EigenMethods.h"
#include "Matrix.h"
#include "Vector.h"

using namespace EigenMethods;

/**
 * parent class for all eigenvalue algorithms
 * @tparam T : vector type
 * @tparam S : matrix type
 */
template <typename T, typename S> class EigenSolver {

private:
    Matrix<S> mMatrix;
    int mMaxIter;
    T mTol;

public :
    /// default constructor
    EigenSolver(const Matrix<S>& matrix, int maxIter, T tol);

    /**
     * @tparam S : matrix type
     * @return the matrix associated to the EigenSolver
     */
    Matrix<S> GetMatrix() ;

    /**
     * @return the tolerance as a double
     */
    double GetTolerance() const;

    /**
     * @return the maximum iteration number as an integer
     */
    int GetMaxIterationNumber() const;
};

#endif //PSC_PROJECT_EIGENVALUES_EIGENSOLVER_H