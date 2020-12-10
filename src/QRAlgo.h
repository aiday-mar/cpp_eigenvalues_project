//
// Created by Aiday on 26/11/2020.
//

#ifndef PSC_PROJECT_EIGENVALUES_POWERALGO_H
#define PSC_PROJECT_EIGENVALUES_POWERALGO_H

#include "EigenSolver.h"
#include "EigenMethods.h"
#include "Matrix.h"
#include "Vector.h"
#include <numeric>
#include <utility>

using namespace EigenMethods;

template<typename T, typename S> class QRAlgo : public EigenSolver<T,S> {
public :

    QRAlgo(Matrix<S> matrix, int maxIter, T tol);

    /**
     * @param A : matrix for which we calculate the QR factorization
     * @return matrices Q and R as a pair, of the QR decomposition
     */
    std::pair<Matrix<S>, Matrix<S>> QRFactorization(Matrix<S> A);

    /**
     * @return a vector of eigenvalues and a matrix containing the corresponding eigenvectors using the matrix stored in the EigenSolver
     */
    std::pair<Vector<T>, Matrix<S>> ComputeEigenvalues();

    /**
     * @param A
     * @return the maximum absolute value in the lower triangular matrix of matrix A
     */
    S maxLowerTriangularMatrix(Matrix<S> A);
};

#endif //PSC_PROJECT_EIGENVALUES_POWERALGO_H
