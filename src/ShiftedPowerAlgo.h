//
// Created by Aiday on 26/11/2020.
//

#ifndef PSC_PROJECT_EIGENVALUES_SHIFTEDPOWERALGO_H
#define PSC_PROJECT_EIGENVALUES_SHIFTEDPOWERALGO_H

#include "EigenSolver.h"
#include "EigenMethods.h"
#include "Matrix.h"
#include "Vector.h"

using namespace EigenMethods;

template<typename T, typename S> class ShiftedPowerAlgo : public EigenSolver<T,S> {
public :
    ShiftedPowerAlgo(const Matrix<S>& matrix, int maxIter, T tol);

    /**
     * @param x_initial
     * @param shift
     * @return the eigenvalue of the matrix closests to the shift value and starting from estimate x_initial
     */
    std::pair<T, Vector<T>> findEigenvalue(const Vector<T>& x_initial, const T shift);
};

#endif //PSC_PROJECT_EIGENVALUES_SHIFTEDPOWERALGO_H
