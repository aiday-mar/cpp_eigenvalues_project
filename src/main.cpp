//
// Created by Moritz on 25/11/20.
//

#include "Matrix.h"

/// main function to use the eigenvalue library
int main(int argc, char* argv[]) {
    // Create a new 2x2 matrix A
    Matrix<double> A(2, 2);

    // Print matrix A
    std::cout << "Print Matrix A:" << std::endl;
    A.Print();

    return 0;
}
