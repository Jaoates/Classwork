# tells user if the inputed string contains the word "accepted"
# None -> None
def read_str():
    s = input("enter a sentance containing 'accepted'")
    while "accepted" not in s:
        s = input("Error: must contain 'accepted'")


# this will compute the factorial of num
# int -> int
def factorial(num):
    num = abs(int(num))
    r = 1
    for i in range(1, num + 1):
        r = r * i
    return r


# this will compute the amount of tax to be paid given an income
# int -> float
def income_tax(income):
    income = int(income)
    if income < 0:
        print("error, income must be positive")
        rate = 0
    elif income < 9951:
        rate = 10
    elif income < 40526:
        rate = 12
    elif income < 86376:
        rate = 22
    elif income < 164926:
        rate = 24
    elif income < 209426:
        rate = 32
    else:
        print("error, income must be less than 209426")
        rate = 0
    return (income * rate) / 100


# takes average of and prints numbers divisible by 7 in interval (inclusive)
# int int -> float
def total_avg(start, end):
    start = int(start)
    end = int(end)
    if end < start:
        print("error, end less than start")
        return 0
    count = 0
    tot = 0
    for i in range(start, end + 1, 1):
        if i % 7 == 0:
            print(i)
            count += 1
            tot += i
    return tot / count
