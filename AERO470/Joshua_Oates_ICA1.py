# 1 Using a while loop, write Python code to compute the factorial of the integer n.
# 2 Use a while loop, with an if statement to print out the even numbers from 0 to 20.  In this activity, use the modulus operator (%)
# 3 Use an if/else statement and the ‘random’ module to code up the rock-paper-scissors game.  Your code should ask the user to pick
# between rock, paper, or scissors, then generate a random choice from the computer.  The user can enter ‘exit’ to quit the game.

from random import *

# 1
def myFactorial(n):
    f = 1
    n += 1
    while n > 1:
        n -= 1
        f*=n
    return f

print(myFactorial(1))
print(myFactorial(2))
print(myFactorial(3))
print(myFactorial(4))
print(myFactorial(0))
print(myFactorial(-1))
print(myFactorial(4.3))

# 2
i = 0
while i <= 20:
    if not i%2:
        print(i)
    i+=1

# 3
def rockPaperScissors():
    print('Rock,Paper,Scissors!')
    print('q to quit')
    while True:
        c = randint(0,3)
        if c == 0:
            c = 'r'
        elif c== 1:
            c = 'p'
        else:
            c = 's'
        u = input("Your Move ('r','p','s'): ")
        if u == 'q':
            return
        elif u not in 'rps':
            print('invalid input')
        elif c == u:
            print(f"computer picked {c}, tie game.")
        elif c == 'r' and u == 'p' or c == 's' and u == 'r' or c ==  'p' and u == 's':
            print(f"computer picked {c}, you win!")
        elif u == 'r' and c == 'p' or u == 's' and c == 'r' or u ==  'p' and c == 's':
            print(f"computer picked {c}, you lose.")

rockPaperScissors()

