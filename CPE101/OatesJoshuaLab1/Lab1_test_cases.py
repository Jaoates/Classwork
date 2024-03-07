"""
 Test cases example for lab 01.

 Name: Joshua Oates

 Section: 07
"""

import unittest  # import the module that supports writing unit testsl


# Define a class that will be used for these test cases.
# This class will include testing functions.
#
# Much of this code should be viewed as boilerplate for now.
#
class TestsLab1(unittest.TestCase):
    def test_expressions(self):  # Defines one testing function.
        self.assertEqual(2 + 3, 5)  # 1
        # Add code after this line (like the line above) for the additional test cases.
        self.assertEqual(2 * 3, 6)  # 2
        self.assertEqual(2 ** 3, 8)  # 3
        self.assertEqual(49 // 3, 16)  # 4
        self.assertAlmostEqual(49 / 3, 16.333333333333332)  # 5
        self.assertEqual(49 / 3.0, 16.333333333333332)  # 6
        self.assertEqual(49 % 3, 1)  # 7
        self.assertEqual(4 * 3 + 17 // 2 - 5, 15)  # 8
        self.assertEqual(4 * (3 + 17) // 2 - 5, 35)  # 9
        self.assertEqual(4 * 3 + (17 // 2) - 5, 15)  # 10
        self.assertEqual(4 * (3 + 17) // (2 - 5), -27)  # 11

    # Run the unit tests.


if __name__ == '__main__':
    unittest.main()

"""
1 How to use interactive mode?
    Type in the prompt at the bottom of the screen.
    
2 Lab1.py - changes applied and run correctly

3 9 test cases for each given expressions

4 Why do we use assertAlmostEqual?
    Floating point numbers may round in such a way that it is not strictly the same
    as you would expect
    
5 difference between / and //
    / is a floating point divide while // is an integer divide. // will only return int
    
6 difference between expression 8, 9, 10, and 11
    although the operators are in the same locations, the parenthesis will change the order of 
    operations, such that they will evaluate to different numbers. 
    In 9, 3+17 will happen "early"
    In 10, 17//2 will happen before 4*3, but this wont change the answer
    In 11, 3+17 and 2-5 will happen "early"

"""
