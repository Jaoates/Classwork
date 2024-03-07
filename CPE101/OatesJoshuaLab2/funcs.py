# Lab 2 Functions
# Name: Joshua Oates
# Section: 07


import math


# 1)purpose statement: This function compute a given formula
# signature: float float -> float


def do_math(x, y):
    res = (x ** 2 - (4 / 3) * x * y ** 2) / ((2 * x ** 2 * y) / 5 - 5)
    return res


# 2)purpose statement: This function will compute the large value of x using inputs
# a, b, and c in a quadratic polynomial 0 = ax^2 + bx + c
# signature: float float float -> float
def quadratic_formula1(a, b, c):
    res = (-b + (b ** 2 - 4 * a * c) ** .5) / (2 * a)
    return res


# 3)purpose statement:This function will compute the small value of x using inputs
# a, b, and c in a quadratic polynomial 0 = ax^2 + bx + c
# signature: float float float -> float
def quadratic_formula2(a, b, c):
    res = (-b - (b ** 2 - 4 * a * c) ** .5) / (2 * a)
    return res


# 4)purpose statement: This function will compute the distance between P1 and P2
# signature: float float float float -> float
def distance(x1, y1, x2, y2):
    res = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
    return res


# 5)purpose statement: returns true if num is negative
# signature: float -> bool
def is_negative(num):
    res = (num < 0)
    return res


# 6)purpose statement: will return true if num is divisible by 5
# signature: float -> bool
def dividable_by_5(num):
    res = 0 == (num % 5)
    return res
