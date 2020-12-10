//
// Created by Moritz on 27/11/20.
//

#include <gtest/gtest.h>
#include "../src/Exception.h"

/// Test if throwing a costume exception works
TEST(testException, testConstructorAndGetters) {
    EXPECT_THROW(throw NegativeNumberException(), NegativeNumberException);
}

