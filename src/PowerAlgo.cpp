//
// Created by Aiday on 26/11/2020.
//

#include "PowerAlgo.h"

template<typename T, typename S>
PowerAlgo<T,S>::PowerAlgo(const Matrix<S>& matrix, int maxIter, T tol) : EigenSolver<T,S>(matrix, maxIter, tol) {
}

template<typename T, typename S>
std::pair<T, Vector<T>> PowerAlgo<T,S>::findLargest(const Vector<T>& x_initial) {

    // inverse of the norm of the initial vector
    T initial_norm_inversed = 1 / x_initial.Norm();
    // normalize x0
    auto x0(x_initial * initial_norm_inversed);

    auto Ax0 = PowerAlgo<T, S>::GetMatrix() * x0;
    double Ax0_norm_inversed = 1 / Ax0.Norm();
    auto x_new(Ax0 * Ax0_norm_inversed);
    // calculating the rayleigh quotient
    T estimated_new_eigenvalue = RayleighQuotient(PowerAlgo<T, S>::GetMatrix(), x_new);

    T estimated_old_eigenvalue = INFINITY;
    auto x_old(x_new * INFINITY);
    int number_iterations = 1;

    // While the difference in the consecutive eigenvalues is bigger than the tolerance
    // and the iteration number is smaller than the maximum iteration number
    while (std::abs(estimated_new_eigenvalue - estimated_old_eigenvalue) >
            PowerAlgo<T, S>::GetTolerance() &&
            number_iterations < PowerAlgo<T, S>::GetMaxIterationNumber()) {

        estimated_old_eigenvalue = estimated_new_eigenvalue;
        x_old = x_new;
        number_iterations += 1;

        auto Ax(PowerAlgo<T, S>::GetMatrix() * x_new);
        double Ax_norm_inversed = 1 / Ax0.Norm();
        x_new = Ax * Ax_norm_inversed;
        estimated_new_eigenvalue = RayleighQuotient(PowerAlgo<T, S>::GetMatrix(), x_new);
    }

    // final eigenvector is normalized
    x_new.Normalize();

    // print message if algo doesn't converge
    if (number_iterations == PowerAlgo<T, S>::GetMaxIterationNumber() &&
           std::abs(estimated_new_eigenvalue - estimated_old_eigenvalue) >
                   PowerAlgo<T, S>::GetTolerance() * estimated_new_eigenvalue) {
        std::cout << "The power method did not converge with " << PowerAlgo<T, S>::GetMaxIterationNumber()
        << " iterations." << std::endl;
    }

    // store in a pair
    return std::make_pair(estimated_new_eigenvalue, x_new);
}

// supporting types (double, float)
template class PowerAlgo<double, double>;
template class PowerAlgo<float, float>;