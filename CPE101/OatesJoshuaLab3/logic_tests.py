# This is a header block example for lab 03.

# You will need to supply the following information.

# Name: Joshua Oates
# Section: 07

import unittest
from logic import *
from conditional import *


# You can dlete pass after wrinting your code
class TestCases(unittest.TestCase):
    def test_is_even(self):
        self.assertTrue(is_even(2))
        self.assertFalse(is_even(3))
        self.assertFalse(is_even(3.3))

    def test_in_an_interval(self):
        # Write test case for the following values: -8, -3, 4, 32, 52, 147.
        self.assertTrue(in_an_interval(-8))
        self.assertFalse(in_an_interval(-3))
        self.assertTrue(in_an_interval(4))
        self.assertFalse(in_an_interval(32))
        self.assertTrue(in_an_interval(52))
        self.assertTrue(in_an_interval(147))

    def test_is_divisible_in_interval(self):
        self.assertTrue(is_divisible_in_interval(70, 7))
        self.assertTrue(is_divisible_in_interval(100, 20))
        self.assertFalse(is_divisible_in_interval(85, 35))
        self.assertFalse(is_divisible_in_interval(40, 20))


    # Run the unit tests.
if __name__ == '__main__':
    unittest.main()
