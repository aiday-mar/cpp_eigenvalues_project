//
// Created by Moritz on 27/11/20.
//

#include "Exception.h"

// Defining constructor for each costume exception
WrongDimensionException::WrongDimensionException() : runtime_error("Dimensions of vectors and/or matrices don't match.") {}
NegativeNumberException::NegativeNumberException() : runtime_error("Input can't be a negative number.") {}
AlreadyAllocatedException::AlreadyAllocatedException() :runtime_error("Memory has already been allocated. Check allocated rows and columns with GetNumberOfRows / GetNumberOfCols.") {}
NonSquareMatrixException::NonSquareMatrixException() : runtime_error("Dimensions of row and columns have to be equal.") {}
SingularMatrixException::SingularMatrixException() : runtime_error("Matrix is singular.") {}
FileNotOpenException::FileNotOpenException() : runtime_error("File is not open") {}