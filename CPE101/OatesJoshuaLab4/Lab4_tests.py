# This is a header block example for lab 03.

# You will need to supply the following information.

# Name:
# Section:

import unittest
from functions import *


# You can dlete pass after wrinting your code
class TestCases(unittest.TestCase):
    def test_factorial(self):
        self.assertEqual(factorial(5), 120)
        self.assertEqual(factorial(.1), 1)
        self.assertEqual(factorial(-4), 24)

    def test_income_tax(self):
        self.assertAlmostEqual(income_tax(9950), 995)
        self.assertAlmostEqual(income_tax(50000), 11000.0)
        self.assertAlmostEqual(income_tax(9951), 1194.12)
        self.assertAlmostEqual(income_tax(87000), 20880.0)
        self.assertAlmostEqual(income_tax(200000), 64000.0)

    def test_total_avg(self):
        self.assertAlmostEqual(total_avg(1,50),28)
        self.assertAlmostEqual(total_avg(50,1),0)
        self.assertAlmostEqual(total_avg(14, 21), 17.5)


# Run the unit tests.
if __name__ == '__main__':
    unittest.main()
