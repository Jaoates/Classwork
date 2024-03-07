# This is a header block example for lab 03.

# You will need to supply the following information.

# Name: Joshua Oates
# Section: 07

# return larger of two numbers
# float float -> float
def max_of_2(num1, num2):
    if num1 > num2:
        return num1
    else:
        return num2



# return largest of 4 numbers
# float float float float -> float
def max_of_4(num1, num2, num3, num4):
    # calls max_of_2
    num1 = max_of_2(num1, num2)
    num2 = max_of_2(num3, num4)
    return max_of_2(num1, num2)


# >= 93   A
# >= 90   A-
# >= 87   B+
# >= 83   B
# >= 80   B-
# >= 77   C+
# >= 73   C
# >= 70   C-
# >= 65   D
# < 65    F
# float -> str

def letter_grade(num):
    if num >= 93:
        out = "A"
    elif num >= 90:
        out = "A-"
    elif num >= 87:
        out = "B+"
    elif num >= 83:
        out = "B"
    elif num >= 80:
        out = "B-"
    elif num >= 77:
        out = "C+"
    elif num >= 73:
        out = "C"
    elif num >= 70:
        out = "C-"
    elif num >= 65:
        out = "D"
    else:
        out = "F"
    return out

