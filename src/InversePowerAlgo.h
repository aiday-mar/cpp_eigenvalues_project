//
// Created by Aiday on 26/11/2020.
//

#ifndef PSC_PROJECT_EIGENVALUES_INVERSEPOWERALGO_H
#define PSC_PROJECT_EIGENVALUES_INVERSEPOWERALGO_H

#include "ShiftedPowerAlgo.h"
#include "EigenMethods.h"
#include "Matrix.h"
#include "Vector.h"

using namespace EigenMethods;

template<typename T, typename S>
class InversePowerAlgo : public ShiftedPowerAlgo<T,S> {
public :
    /// default constructor
    InversePowerAlgo(const Matrix<S>& matrix, int maxIter, T tol);

    /**
     * @param x_initial
     * @return the smallest eigenvalue of the matrix and the corresponding eigenvector where you use the initial estimate x_initial
     */
    std::pair<T, Vector<T>> findSmallest(const Vector<T>& x_initial);
};

#endif //PSC_PROJECT_EIGENVALUES_INVERSEPOWERALGO_H
