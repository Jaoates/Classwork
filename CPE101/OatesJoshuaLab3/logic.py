# This is a header block example for lab 03.

# You will need to supply the following information.

# Name: Joshua Oates
# Section: 07

# tells if an int is even
# int -> bool
def is_even(num):
    return not bool(num % 2)


#  [-12,-6), [-2, 5], (15,28), (42, 52], and [110,237].
# float -> bool
def in_an_interval(num):
    i1 = bool(-12 <= num and -6 > num)
    i2 = bool(-2 <= num and 5 >= num)
    i3 = bool(15 < num and 28 > num)
    i4 = bool(42 < num and 52 >= num)
    i5 = bool(110 <= num and 237 >= num)
    return bool(i1 or i2 or i3 or i4 or i5)


# returns True when the 1st argument is in the (60, 140] interval and is divisible to 2nd argument, which is in the interval of [-10,40].
# int int -> bool
def is_divisible_in_interval(num1, num2):
    i1 = bool(60 < num1 and 140 >= num1)
    i2 = bool(-10 <= num2 and 40 >= num2)
    i3 = not bool(num1 % num2)
    return bool(i1 and i2 and i3)
