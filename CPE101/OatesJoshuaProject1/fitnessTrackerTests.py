# Project 1 Fitness Tracking Tests
# At least 2 tests for each function
# Name: Joshua Oates
# Section: 07
# Date: 4/30

import unittest
from fitnessTrackerFuncs import *


# compute_BMI()
# compute_MET()
# compute_equivalent_miles()
# compute_calorie_burned()
# BMI_category()

class MyTestCase(unittest.TestCase):
    def test_convert_lb_to_kg(self):
        self.assertAlmostEqual(convert_lb_to_kg(2), 0.90718474)
        self.assertAlmostEqual(convert_lb_to_kg(0), 0)
        self.assertAlmostEqual(convert_lb_to_kg(178), 80.73944186)

    def test_compute_BMI(self):
        self.assertAlmostEqual(compute_BMI(140, 72), 18.98533950617284)
        self.assertAlmostEqual(compute_BMI(179, 72), 24.274112654320987)
        self.assertAlmostEqual(compute_BMI(140, 60), 27.33888888888889)

    def test_compute_MET(self):
        self.assertEqual(compute_MET(1), 5)
        self.assertEqual(compute_MET(2), 7)
        self.assertEqual(compute_MET(3), 3)
        self.assertEqual(compute_MET(4), 4)

    def test_equivalent_miles(self):
        self.assertEqual(compute_equivalent_miles(72, 1, 120), 6.758181818181817)
        self.assertEqual(compute_equivalent_miles(72, 3, 120), 5.631818181818181)
        self.assertEqual(compute_equivalent_miles(60, 3, 50), 1.955492424242424)

    def test_calorie_burned(self):
        self.assertEqual(compute_calorie_burned(20, 80, 7), 196.0)
        self.assertEqual(compute_calorie_burned(60, 80, 7), 588.0)
        self.assertEqual(compute_calorie_burned(15, 110, 5), 144.375)

    def test_BMI_category(self):
        self.assertEqual(BMI_category(18), 'Underweight')
        self.assertEqual(BMI_category(28), 'Over Weight')
        self.assertEqual(BMI_category(66), 'Obesity')


if __name__ == '__main__':
    unittest.main()
