# PSCS Project, 2020, Eigenvalues

## Authors

- Aiday Marlen Kyzy
- Carmine Schipani
- Moritz Waldleben

## Abstract

In this project, we created a simple c++ library to calculate eigenvalues and eigenvectors. The library is build with a CMakeLists.txt created in CLion and uses googletests for the testing. A templated Matrix and Vector class provides functions and operator overloading such that matrices and vectors can be used in an natural way. Supported are int, float, double types. The library overs different kinds of so called EigenSolvers such as the power method, the inverse power method and the shifted inverse power method. Furthermore a QR methods is implemented as well. A documentation is provided in the docs with doxygen.

## Library structure
The structure is shown below. For simplicity the corresponding hpp files are omited. They are placed in the same folder as the corresponding cpp files.

    pcsc_project_eigenvalues
    ├── googletest
    ├── src
    ├── ├── main.cpp
    ├── ├── Exceptions.cpp
    ├── ├── Matrix.cpp
    ├── ├── Vector.cpp
    ├── ├── EigenMethods.cpp
    ├── ├── EigenSolver.cpp
    ├── ├── PowerAlgo.cpp
    ├── ├── InversPowerAlgo.cpp
    ├── ├── ShiftePowerAlgo.cpp
    ├── ├── QRAlgo.cpp
    ├── test
    ├── ├── testException.cpp
    ├── ├── testMatrix.cpp
    ├── ├── testVectors.cpp
    ├── ├── testEigenMethods.cpp
    ├── ├── testAlgos.cpp
    ├── ├── QRAlgo.cpp
    ├── CMakeList.txt
    ├── README.md


## Use the Library

The CMakeList provides two target files:

`run_eigenvalues`
To use the library use the main.cpp file and do the calculations in the main function. A basic example is provided there. All source files are grouped in a library eigenvalues which could than be used in external project.

`test_eigenvalus`
Testing files are build here. In CLion the googletest framework can be accesed directly and run seperatly.

To see the program execution you may run for example testAlgos.cpp. This test creates an instance of an EigenSolver subclass such as PowerAlgo, and runs its' method findLargest to find the largest eigenvalue which is then printed in the console. In order to test the library, you may add your own example to the main.cpp file and run the project. Below you may find two such examples. The first example creates manually a Matrix, and the second imports it from a text file. 

## Two simple examples
1. Create a simple matrix and calculate the smallest eigenvalues
```
// Initialize a 2x2 identity matrix
Matrix<double> A(2, 2);
M(0, 0) = 0, Matrix(0, 1) = 1;
M(1, 0) = 1, Matrix(1, 1) = 0;

// Run the power method
int maxiter = 100;
double tol = 0,01;
Vector<double> x0(2); // x0 = [0, 0] by default

PowerAlgo<double, double> p_algo = PowerAlgo<double, double>(A, maxiter, tol);
std::pair<double, Vector<double>> eigenval_and_vect = p_algo.findLargest(x0);

// Output largest eigenvalue
std::cout << std::get<0>(eigenval_and_vect)
```

2. Load matrix from specified file path and calculate the QR factorization
```
// Load matrix form file path
Matrix<double> B();
B.LoadFromTXTFile("home/user/downloads/B.txt")

// Run QR algorithm
int max_iter = 100;
double tol = 0.01;
Vector<double> x0(2); // x0 = [0, 0] by default

QRAlgo<double, double> qr_algo = QRAlgo<double, double>(B, max_iter, tol);
std::pair<Matrix<double>, Matrix<double>> qr = qr_algorithm_instance.QRFactorization(A);

// Print Q
std::cout << "Q:" << std::endl << get<0>(qr);
```

## List of features

The classes we implemented can already be seen in the library structure diagram above. We explain below for the sake of completeness the classes we have implemented. 

**Exceptions.cpp**

A class implementing exceptions specific to our project. It contains exceptions such as `WrongDimensions`.

**Vector.cpp**

A class implementing vectors. Initially we used the standard library implementation `std::vector<T>` but decided to implement this class in order to be coherent with the Matrix class we defined previously. This vector class implements operators and vector specific operations such as finding the norm of a vector. 

**Matrix.cpp**

The first class we have implemented. This class is similar in its construction to the vector class, it implements operators and methods such as taking the inverse, finding the determinant, to cite a few examples. 

**EigenMethods.cpp**

This file contains methods which are used frequently throughout the project but have not been placed into a class. It includes a method to compute the Rayleigh quotient, and one to compute the projection of a vector on another vector. 

**EigenSolver.cpp**

This is the base class from which our main classes to calculate the eigenvalues of a matric inherit from. It is not instantiated as such in the tests. It contains as parameters the matrix for which we want to find the eigenvalues, the maximum number of iterations and the tolerance used in the subclasses. 

**PowerAlgo.cpp**

This class implements the power method to find the largest eigenvalue of a matrix.

**InversePowerAlgo.cpp**

This class implements the inverse power method to find the smallest eigenvalue of a matrix. It runs essentially the ShiftedPowerAlgo method by setting the shift equal to zero. 

**ShiftedPowerAlgo.cpp**

This class implements the shifted power algo method. It finds the eigenvalue of a matrix closest to a certain shift value. 

**QRAlgo.cpp**

This class implements the QR algorithm to find all the eigenvalues. It contains a method to calculate the QR decomposition of a matrix. This method is later used to return the eigenvalues in a Vector, and the corresponding eigenvectors in a Matrix in the findEigenvalues method. 

## Possible Improvements

First it should be noted that the methods we have implemented return the eigenvalues only in the case when these are real. We have not had the time to implement these same methods with complex eigenvalues. However we have tested the Matrix.cpp and Vector.cpp files with std::complex<T> and know that their methods work with the complex type. Next we have previously run into problems with making the Vector inherit from the Matrix, and it was in the end decided that these would be defined as separate classes. This means that we could have written less code for the Vector class. 

Aditionally, while we have tried our best to follow the guidelines for concise structured code implementing polymorphism in c++, there may be certain classes where the code could have been structured better. Finally in order to facilitate user interaction with the library, we could create a GUI interface, in order to have to avoid dealing with the console.

## Reflection on the Project

The project has been a pleasure to work on. We have learned a lot, gained valuable skills in c++ and learned to work as a team, by splitting the work. Aside from the knowledge in c++, we have also learned how to develop a small project from beginning to the end, and have been familiarized with new technologies such as GitLab, Doxygen and Google Test Suites. The hard part was most likely the debugging, and spending hours finding small errors that had propagated in the whole project. However with hard and meticulous work, we have been able to finally deliver the project. 
