//
// Created by Moritz on 27/11/20.
//

#ifndef PSC_PROJECT_EIGENVALUES_EXCEPTION_H
#define PSC_PROJECT_EIGENVALUES_EXCEPTION_H

#include <string>
#include <stdexcept>

/// Costume Exceptions
struct WrongDimensionException : public std::runtime_error {
    WrongDimensionException();
};

/// Costume Exceptions
struct NegativeNumberException : public std::runtime_error{
    NegativeNumberException();
};

/// Costume Exceptions
struct AlreadyAllocatedException : public std::runtime_error {
    AlreadyAllocatedException();
};

/// Costume Exceptions
struct NonSquareMatrixException : public std::runtime_error{
    NonSquareMatrixException();
};

/// Costume Exceptions
struct SingularMatrixException : public std::runtime_error{
    SingularMatrixException();
};

/// Costume Exceptions
struct FileNotOpenException : public std::runtime_error {
    FileNotOpenException();
};

#endif //PSC_PROJECT_EIGENVALUES_EXCEPTION_H

