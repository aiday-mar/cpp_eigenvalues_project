//
// Created by Aiday on 23/11/2020.
//

#ifndef PSC_PROJECT_EIGENVALUES_EIGENMETHODS_H
#define PSC_PROJECT_EIGENVALUES_EIGENMETHODS_H

#include "Matrix.h"
#include "Vector.h"
#include <complex>
#include <cmath>

/// a namespace which groups methods which are used in the EigenSolver
namespace EigenMethods {

    /**
     * @param A : matrix used in the Rayleigh quotient
     * @param x : vector used in the Rayleigh quotient
     * @return the rayleigh quotient with A and x
     */
    template<typename T, typename S>
    T RayleighQuotient(Matrix<S> A, Vector<T> x);

    /**
     * @param a : vector a which is projected on u
     * @param u : vector on which projection is done
     * @return the projection of a on u
     */
    template<typename T>
    Vector<T> projection(Vector<T> u, Vector<T> a);

}

#endif //PSC_PROJECT_EIGENVALUES_EIGENMETHODS_H
