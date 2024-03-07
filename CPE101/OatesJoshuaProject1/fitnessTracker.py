# Main Program
# Project 1 Fitness Tracking
#
# Name: Joshua Oates
# Section: 07
# Date: 4/30
from fitnessTrackerFuncs import *


def main():
    welcome()
    another_user = "y"
    while another_user == "y":
        name = input_name()
        age = input_age()
        height = input_height()
        weight = input_weight()
        duration = input_duration()
        e_type = input_exercise_type()
        met_value = compute_MET(e_type)
        miles = compute_equivalent_miles(height, e_type, duration)
        calories = compute_calorie_burned(duration, convert_lb_to_kg(weight), met_value)
        BMI = compute_BMI(weight, height)
        BMIC = BMI_category(BMI)
        print_info(name, age, height, weight, calories, BMI, miles, BMIC)
        another_user = input("Do you want to apply fitness tracking for another user (y/n)? ")


if __name__ == '__main__':
    main()
