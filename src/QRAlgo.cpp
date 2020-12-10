//
// Created by Aiday on 26/11/2020.
//

#include "QRAlgo.h"

template<typename T, typename S>
QRAlgo<T,S>::QRAlgo(Matrix<S> matrix, int maxIter, T tol) : EigenSolver<T,S>(matrix, maxIter, tol) {}


template<typename T, typename S>
std::pair<Matrix<S>, Matrix<S>> QRAlgo<T,S>::QRFactorization(Matrix<S> A){

    // initialized matrices Q, R and U
    // where U will contain intermediate vectors used in the calculation of the QR method
    Matrix<S> Q(A.GetNumberOfRows(), A.GetNumberOfCols());
    Matrix<S> R(A.GetNumberOfRows(), A.GetNumberOfCols());
    Matrix<S> U(A.GetNumberOfRows(), A.GetNumberOfCols());

    for(int i= 0; i<A.GetNumberOfRows(); i++) {
        for(int j= 0; j<A.GetNumberOfCols(); j++){
            R(i,j) = 0;
        }
    }

    // the 1st column of matrix U is the 1st column of matrix A
    Vector<S> old_column = A.getColumn(0);
    for(int i=0; i< A.GetNumberOfRows(); i++) {
        U(i,0) = old_column(i);
    }

    T norm0 = 1/old_column.Norm();
    auto column0(old_column * norm0);

    // the 1st column of matrix Q is the normalized 1st column of matrix A
    for(int i=0; i< A.GetNumberOfRows(); i++) {
        Q(i,0) = column0(i);
    }

    Vector<S> a = Q.getColumn(0);
    Vector<S> b = A.getColumn(0);
    R(0,0) = a*b;

    for(int j=1; j< A.GetNumberOfCols(); j++) {

        Vector<S> column = A.getColumn(j);

        for(int k=0; k < j;k++){
            Vector<S> columnU = U.getColumn(k);
            // jth column of A is projected onto columnU
            Vector<S> projectionVec = projection(column,columnU);
            // remove that column from the jth column of A
            column = column - projectionVec;
        }

        // intermediate vector stored in the jth column of U
        for(int i=0; i< A.GetNumberOfRows(); i++) {
            U(i,j) = column(i);
        }

        // after normalized it can be added to matrix Q
        column = column*(1/column.Norm());

        for(int i=0; i< A.GetNumberOfRows(); i++) {
            Q(i,j) = column(i);
        }

        // updating matrix Q
        // <= because we include the diagonal
        for(int i=0; i<=j;i++) {
            Vector<S> a = Q.getColumn(i);
            Vector<S> b = A.getColumn(j);
            R(i,j) = a*b;
        }
    }

    // return matrices Q and R in a pair
    std::pair<Matrix<S>, Matrix<S>> pair = std::make_pair(Q,R);
    return pair;
}

template <typename T, typename S>
std::pair<Vector<T>, Matrix<S>>  QRAlgo<T,S>::ComputeEigenvalues(){

    int number_of_iterations = 0;
    std::pair<Matrix<S>, Matrix<S>> pair = QRFactorization(QRAlgo<T, S>::GetMatrix());
    // multiply the matrix 1 and 0 in the inverse order
    Matrix<S> A = std::get<1>(pair)*std::get<0>(pair);
    S maxLower = maxLowerTriangularMatrix(A);

    // we require the maximum in the lower triangle to be lower than a specific tolerance
    while((number_of_iterations < QRAlgo<T, S>::GetMaxIterationNumber()) && (maxLower > QRAlgo<T, S>::GetTolerance())){
        // continue to calculate the QR pair
        pair = QRFactorization(A);
        // continue multiplying the Q and R in inverse order
        A = std::get<1>(pair)*std::get<0>(pair);
        maxLower = maxLowerTriangularMatrix(A);
        number_of_iterations += 1;
    }

    Vector<S> approximateEigenvalues(QRAlgo<T, S>::GetMatrix().GetNumberOfRows());
    // assume the matrix is diagonal
    for(int i = 0; i < QRAlgo<T, S>::GetMatrix().GetNumberOfRows(); i++){
        approximateEigenvalues(i) = A(i,i);
    }

    std::pair<Vector<S>, Matrix<S>> finalPair = std::make_pair(approximateEigenvalues,A);
    return finalPair;
}

template<typename T, typename S>
S QRAlgo<T,S>::maxLowerTriangularMatrix(Matrix<S> A){

    // current max is zero
    T max = 0;
    for(int i=0; i < A.GetNumberOfRows(); i++) {
        for(int j=0; j < i; j++) {
            // suppose the absolute value of an element in the strict lower triangle is bigger than max
            if(std::abs(A(i,j)) > max){
                // then this current max is set to this value
                max = std::abs(A(i,j));
            }
        }
    }
    return max;
}

// supported types (double, float)
template class QRAlgo<double, double>;
template class QRAlgo<float, float>;
