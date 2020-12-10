cmake_minimum_required(VERSION 3.17)
project(psc_project_eigenvalues)

set(CMAKE_CXX_STANDARD 14)

add_subdirectory(googletest)

add_library(eigenvalues src/Exception.cpp src/Exception.h src/Matrix.cpp src/Matrix.h src/Vector.cpp src/Vector.h
        src/EigenMethods.cpp src/EigenMethods.h src/EigenSolver.cpp src/EigenSolver.h src/PowerAlgo.cpp src/PowerAlgo.h
        src/ShiftedPowerAlgo.cpp src/ShiftedPowerAlgo.h src/InversePowerAlgo.cpp src/InversePowerAlgo.h src/QRAlgo.cpp src/QRAlgo.h)

add_executable(run_eigenvalues src/main.cpp)
target_link_libraries(run_eigenvalues eigenvalues)

add_executable(test_eigenvalues test/testException.cpp test/testMatrix.cpp test/testVector.cpp test/testEigenMethods.cpp test/testAlgos.cpp test/testQRAlgo.cpp)
target_link_libraries(test_eigenvalues gtest_main gtest pthread eigenvalues)