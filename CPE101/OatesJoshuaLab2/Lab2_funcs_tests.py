# Lab 2 Test Cases
# Name: Joshua Oates
# Section: 07
##############################################################
import unittest
import funcs


# 3 test cases for each function
class TestCases(unittest.TestCase):
    def test_do_math(self):
        self.assertAlmostEqual(funcs.do_math(2.3, 6), -13.65774428)
        self.assertAlmostEqual(funcs.do_math(1, 0), -0.2)
        self.assertAlmostEqual(funcs.do_math(-66, 1000), 50.50769544)

    def test_quadratic_formula1(self):
        self.assertAlmostEqual(funcs.quadratic_formula1(2, -10, 4), 4.561552813)
        self.assertAlmostEqual(funcs.quadratic_formula1(-1, 4, -1), 0.2679491924)
        self.assertAlmostEqual(funcs.quadratic_formula1(-1, 0, 0), 0)

    def test_quadratic_formula2(self):
        self.assertAlmostEqual(funcs.quadratic_formula2(2, -10, 4), .4384471872)
        self.assertAlmostEqual(funcs.quadratic_formula2(-1, 4, -1), 3.732050808)
        self.assertAlmostEqual(funcs.quadratic_formula2(-1, 0, 0), 0)

    def test_distance(self):
        self.assertAlmostEqual(funcs.distance(1, 0, 3, 0), 2)
        self.assertAlmostEqual(funcs.distance(1, 1, 3, 3), 2.8284271247461903)
        self.assertAlmostEqual(funcs.distance(10, -13131, 1, 999), 14130.002866241748)

    def test_is_negative(self):
        self.assertFalse(funcs.is_negative(5))
        self.assertTrue(funcs.is_negative(-5))
        self.assertFalse(funcs.is_negative(0))

    def test_dividable_by_5(self):
        self.assertTrue(funcs.dividable_by_5(5))
        self.assertTrue(funcs.dividable_by_5(20000000))
        self.assertFalse(funcs.dividable_by_5(5053661))


# Run the unit tests.
if __name__ == '__main__':
    unittest.main()
