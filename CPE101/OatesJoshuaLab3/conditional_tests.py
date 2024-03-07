# This is a header block example for lab 03.

# You will need to supply the following information.

# Name: Joshua Oates
# Section: 07

import unittest
from conditional import *


# You can delete pass after wrinting your code
class TestCases(unittest.TestCase):
    def test_max_of_2(self):
        self.assertEqual(max_of_2(3,4),4)
        self.assertEqual(max_of_2(-2,-4),-2)
        self.assertEqual(max_of_2(3,6),6)
    def test_max_of_4(self):
        self.assertEqual(max_of_4(1,2,3,4),4)
        self.assertEqual(max_of_4(1,2,3,5),5)
        self.assertEqual(max_of_4(-1,-2,-3,-1),-1)
    def test_letter_name(self):
        self.assertEqual(letter_grade(91),"A-")
        self.assertEqual(letter_grade(89.999),"B+")
        self.assertEqual(letter_grade(94),"A")
        self.assertEqual(letter_grade(70),"C-")
        self.assertEqual(letter_grade(72),"C-")
        self.assertEqual(letter_grade(66),"D")
        self.assertEqual(letter_grade(65),"D")
        self.assertEqual(letter_grade(64),"F")
        self.assertEqual(letter_grade(10),"F")





# Run the unit tests.
if __name__ == '__main__':
    unittest.main()
