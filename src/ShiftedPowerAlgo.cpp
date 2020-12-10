//
// Created by Aiday on 26/11/2020.
//

#include "ShiftedPowerAlgo.h"

template<typename T, typename S>
ShiftedPowerAlgo<T,S>::ShiftedPowerAlgo(const Matrix<S>& matrix, int maxIter, T tol) : EigenSolver<T,S>(matrix, maxIter, tol) {
}

template<typename T, typename S>
std::pair<T, Vector<T>> ShiftedPowerAlgo<T,S>::findEigenvalue(const Vector<T>& x_initial, const T shift) {

    // inverse of the norm of the initial vector
    T initial_norm_inversed = 1 / x_initial.Norm();
    // normalize x0
    auto x0(x_initial * initial_norm_inversed);

    // diagonal matrix with shift on the diagonal
    Matrix<T> ShiftId(x0.GetSize(), x0.GetSize());
    for(int i = 0; i < x0.GetSize(); i++) {
        ShiftId(i, i) = shift;
    };

    // remove the shifted matrix from the matrix and multiply by the normalized initial vector
    auto Ax0 = (ShiftedPowerAlgo<T, S>::GetMatrix() - ShiftId).Inverse() * x0;
    T Ax0_norm_inversed = 1 / Ax0.Norm();
    auto x_new(Ax0 * Ax0_norm_inversed);
    // compute the corresponding rayleigh quotient
    T estimated_new_eigenvalue = RayleighQuotient(ShiftedPowerAlgo<T, S>::GetMatrix(), x_new);

    T estimated_old_eigenvalue = INFINITY;
    auto x_old(x_new * INFINITY);
    int number_iterations = 1;

    // while the difference in the consecutive eigenvalues is bigger than the tolerance
    // and the iteration number is smaller than the maximum iteration number
    while (std::abs(estimated_new_eigenvalue - estimated_old_eigenvalue) >
                   ShiftedPowerAlgo<T, S>::GetTolerance() &&
           number_iterations < ShiftedPowerAlgo<T, S>::GetMaxIterationNumber()) {

        estimated_old_eigenvalue = estimated_new_eigenvalue;
        x_old = x_new;
        number_iterations += 1;

        auto Ax((ShiftedPowerAlgo<T, S>::GetMatrix() - ShiftId).Inverse() * x_new);
        T Ax_norm_inversed = 1 / Ax0.Norm();
        x_new = Ax * Ax_norm_inversed;
        estimated_new_eigenvalue = RayleighQuotient(ShiftedPowerAlgo<T, S>::GetMatrix(), x_new);
    }

    // normalize the eigenvector
    x_new.Normalize();

    // print message if algo doesn't converge
    if (number_iterations == ShiftedPowerAlgo<T, S>::GetMaxIterationNumber() &&
        std::abs(estimated_new_eigenvalue - estimated_old_eigenvalue) >
        ShiftedPowerAlgo<T, S>::GetTolerance() * estimated_new_eigenvalue) {
        std::cout << "The inverse or shifted inverse method did not converge with " << ShiftedPowerAlgo<T, S>::GetMaxIterationNumber()
                  << " iterations." << std::endl;
    }

    // return the eigenvalue and the corresponding eigenvector as a pair
    return std::make_pair(estimated_new_eigenvalue, x_new);
}

// supported types (double, float)
template class ShiftedPowerAlgo<double, double>;
template class ShiftedPowerAlgo<float, float>;