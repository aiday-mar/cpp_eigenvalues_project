//
// Created by Aiday on 26/11/2020.
//

#include "InversePowerAlgo.h"

template<typename T, typename S>
InversePowerAlgo<T,S>::InversePowerAlgo(const Matrix<S> &matrix, int maxIter, T tol) : ShiftedPowerAlgo<T,S>(matrix, maxIter, tol) {
}

template<typename T, typename S>
std::pair<T, Vector<T>> InversePowerAlgo<T,S>::findSmallest(const Vector<T> &x_initial) {
    // with zero shift it's the inverse method
    T shift = 0;
    return ShiftedPowerAlgo<T, S>::findEigenvalue(x_initial, shift);
}

// supported types (double, float)
template class InversePowerAlgo<double, double>;
template class InversePowerAlgo<float, float>;