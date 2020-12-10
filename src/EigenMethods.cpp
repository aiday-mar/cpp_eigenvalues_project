//
// Created by Aiday on 23/11/2020.
//

#include "EigenMethods.h"

namespace EigenMethods {

    template<typename T, typename S>
    T RayleighQuotient(Matrix<S> A, Vector<T> x) {

        Vector<S> temp(A * x);
        auto size = temp.GetSize();
        T numerator = 0;
        T denominator = 0;

        // performing calculations to find the numerator and the denominator values
        for (int i = 0; i < size; i++) {
            numerator = numerator + x(i) * temp(i);
            denominator = denominator + x(i) * x(i);
        }

        return numerator / denominator;
    }

    template<typename T>
    Vector<T> projection(Vector<T> a, Vector<T> u) {

        // scalar product of u and a
        double numerator = u * a;
        // scalar product of u and u
        double denominator = u * u;
        double fraction = numerator / denominator;

        // vector res is initialized to be the projection of a on u
        Vector<T> res(u * fraction);
        return res;
    }

    // supported types (double, float)
    template double RayleighQuotient(Matrix<double> A, Vector<double> x);
    template float RayleighQuotient(Matrix<float> A, Vector<float> x);

    template Vector<double> projection(Vector<double> u, Vector<double> a);
    template Vector<float> projection(Vector<float> u, Vector<float> a);
}