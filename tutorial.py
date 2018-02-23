###Python works as a calculator
1+1

7-4

4*2

4/2

5/2

5%2

5./2

float(5)/2

###Data types

print(Hello world!)

print("Hello world!")

"5."/2

type("5.")

type(5)

type(5.)

###Defining variables

x = 1

1 + x

y = "um"

###Boolean Logic
x == y

x > 0

x <= 1

if x > 0:
	print("Do something")

elif x < 0:
	print("Do something else")

if x > 0 and x == 2:
	print("x is 2")

else:
	x = 2

if x != 2:
	print("x is not 2")

if x > 0 and x == 2:
	print("x is 2")

###Data storers
x = 1
y = 2
z = 3

numbers = [1,2,3]

print(x)
print(numbers[0])

numbers_2 = numbers + [x,y,z]
print(numbers)
print(numbers_2)


JEST = ["J","E","S","T"]
ANBIOQ = ["A","N","B","I","O","Q"]

###Loops
for letter in ANBIOQ:
	print(letter)

###User-defined functions

x = 4 + 4

def soma(number_1, number2):
	final_number = number_1 + number_2
	return final_number

y = soma(5,3)

###Libraries
import math

z = sum(5,3)
