# Project 1 Fitness Tracking
#
# Name: Joshua Oates
# Section: 07
# Date: 4/30

# purpose: Display the welcome message.
# signature: None -> None
def welcome():
    """
    display the following:
        'Welcome to Fitness Tracker Application!

            To begin you must specify the participant's name, age, height (in inches),
            weight (in lb), exercise duration (in minutes), and the exercise type (1-4)!
            Now you can compute the burned calories and BMI.'
    """

    print("Welcome to Fitness Tracker Application!")
    print()
    print("    To begin you must specify the participant's name, age, height (in inches),")
    print("    weight (in lb), exercise duration (in minutes), and the exercise type (1-4)!")
    print("    Now you can compute the burned calories and BMI.")
    pass


# purpose: This function has no parameters and returns a string. This function must prompt a user for the name.
# signature: None -> string
def input_name():
    """
    1) prompt user with:
        'Enter the name: '
    2) get the name
    3) return the name
    """
    return input('Enter the name: ')


# purpose: This function has no parameters and returns a positive integer.
# signature: None -> int
def input_age():
    """
    1) prompt user with:
        'Enter age: '
    2) get age
    3) if age > 0
    then
    return age
    4) else
    print 'Error: age must be an integer > 0!'
    goto 1)
    """
    while 1:
        age = int(input('Enter age: '))
        if age > 0:
            return age
        else:
            print('Error: age must be an integer > 0!')


# purpose: This function has no parameters and returns a positive floating point number as height in inches.
# signature: None -> float
def input_height():
    """
    1) prompt user with:
        'Enter height in inches: '
    2) get height
    3) if height > 0
    then
    return height
    4) else
    print 'Error: height must be greater than  0!'
    goto 1)
    """
    while 1:
        height = float(input('Enter height in inches: '))
        if height > 0:
            return height
        else:
            print('Error: height must be greater than  0!')


# purpose:his function has no parameters and return a positive floating point number as weight in pound (lb).
# signature: None -> float
def input_weight():
    """
    1) prompt user with:
        'Enter weight in lb: '
    2) get weight
    3) if weight > 0
    then
    return weight
    4) else
    print 'Error: weight must be greater than 0!'
    goto 1)
    """
    while 1:
        weight = float(input('Enter weight in lb: '))
        if weight > 0:
            return weight
        else:
            print('Error: weight must be greater than 0!')


# purpose: This function has no parameters and return a positive floating point number. This function must prompt a
# user for a positive real value and return it.
# signature: None -> float
def input_duration():
    """
    1) prompt user with:
        'Enter duration of exercise in minutes: '
    2) get duration
    3) if duration > 0
    then
    return weight
    4) else
    print 'Error: duration must be greater than 0!'
    goto 1)
    """
    while 1:
        duration = float(input('Enter duration of exercise in minutes: '))
        if duration > 0:
            return duration
        else:
            print('Error: duration must be greater than 0!')


# purpose: get exercise type
# signature: None -> int
def input_exercise_type():
    """
    1) prompt user with:
        'Enter exercise type: 1) low-impact 2) high-impact 3) slow-paced 4) fast-paced '
    2) get type
    3) if type is an int in [1,4]
    then
    return weight
    4) else
    print 'Error: exercise type must be an integer between [1,4]! '
    goto 2)
    reprompt with:
        'Enter exercise type [1,4]:'
    6) goto 2)
    """
    e_type = int(input('Enter exercise type: 1) low-impact 2) high-impact 3) slow-paced 4) fast-paced '))
    if e_type <= 4 and e_type >= 1:
        return e_type
    else:
        print('Error: exercise type must be an integer between [1,4]!')
    while 1:
        if e_type <= 4 and e_type >= 1:
            return e_type
        else:
            e_type = int(input('Enter exercise type [1,4]: '))


# purpose:print_info(name, age, height, weight, calorie_burned, bmi, miles) – This function displays the
# information in the following order the floating point numbers must have 2 digits after decimal
# point. (Check the format in sample output)
# signature: str, int, float, float, float, float, float, -> None
def print_info(name, age, height, weight, calories, BMI, miles, BMIC):
    """
    print:
               Name :
                Age :
             Height :
             Weight :
              Miles :
    Burned Calories :
                BMI :
       BMI Category :
    """
    print("           Name :", end=" ")
    print(name)
    print("            Age :", end="")
    print("%5d"%age)
    print("         Height :", end="")
    print("%7.2f"%height)
    print("         Weight :", end="")
    print("%7.2f"%weight)
    print("          Miles :", end="")
    print("%7.2f"%miles)
    print("Burned Calories :", end="")
    print("%7.2f"%calories)
    print("            BMI :", end="")
    print("%7.2f"%BMI)
    print("   BMI Category :", end=" ")
    print(BMIC)


# purpose: converts lbs to kg
# signature: float -> float
def convert_lb_to_kg(weight):
    """
    1) get weight
    2) return weight * 0.45359237
    """
    return weight * 0.45359237


# purpose: computes BMI
# signature: int -> int
def compute_MET(exercise_type):
    """
    1) get exercise_type
    2) if exercise_type = 1
    then
        return 5
    3) if exercise_type = 2
    then
        return 7
    4) if exercise_type = 3
    then
        return 3
    5) if exercise_type = 4
    then
        return 4
    """
    if exercise_type == 1:
        return 5
    elif exercise_type == 2:
        return 7
    elif exercise_type == 3:
        return 3
    else:
        return 4


# purpose: Total calories burned = Duration (in minutes)*(MET*3.5*weight in kg)/200
# signature: float float int -> float
def compute_calorie_burned(duration, weight, met_value):
    """
    1)return Total calories burned = Duration (in minutes)*(MET*3.5*weight in kg)/200
    """
    return duration * (met_value * 3.5 * weight) / 200


# purpose: computes BMI from weight in lb and height in in
# signature: float float -> float
def compute_BMI(weight, height):
    """
    1) get weight height
    2) return BMI = 703*weight/height**2
    """
    return 703 * weight / height ** 2


# purpose: gives BMI catagory
# signature: float -> str
def BMI_category(BMI):
    """
    1) get BMI
    2) if BMI < 18.59
        cata = "Underweight"
    3) elif BMI < 25
        cata = "Normal Weight"
    4) elif BMI < 30
        cata = "Over Weight"
    5) else
        cata = "Obesity"
    6) return cata
    """
    if BMI < 18.59:
        c = "Underweight"
    elif BMI < 25:
        c = "Normal Weight"
    elif BMI < 30:
        c = "Over Weight"
    else:
        c = "Obesity"
    return c


# purpose: takes height exercise_type and duration and returns equivalent miles
# signature: float int float -> float
def compute_equivalent_miles(height, exercise_type, duration):
    """
    1) get exercise_type, height, duration
    2) if exercise_type = 1
    then
        steps = 120
    3) if exercise_type = 2
    then
        steps = 227
    4) if exercise_type = 3
    then
        steps = 100
    5) if exercise_type = 4
    then
        steps = 152
    7) milesInStep = ((.413*height)/12)/5280
    8) return milesInStep*steps*duration
    """

    if exercise_type == 1:
        steps = 120
    elif exercise_type == 2:
        steps = 227
    elif exercise_type == 3:
        steps = 100
    else:
        steps = 152
    miles_in_step = ((.413 * height) / 12) / 5280
    return miles_in_step * steps * duration
