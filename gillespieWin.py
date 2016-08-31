from ctypes import *

# This line causes the error (Keagan)
lib = cdll.LoadLibrary('./libfoo.so')

##f.bar(c_double(1.2))
#lib.f.restype = c_float
#print lib.f(c_float(1.3))

charty = c_char_p
charptr = POINTER(charty)
lib.gillespie.restype = charptr
lib.gillespie.argtypes = [charptr]

# math model is defined in c codes. 
# three species: init copy numbers = [10, 0, 100]
# 4 parameters: 
# example: species 1 synthesis rate, deg rate = 0.4, 2
# Third lasts param: time interval
# Second last param: Number of data points.
# the last param: not necessary at this point.
param_input = c_char_p("10 0  100  0.4    2   100   0.1  0.02   10   240000")
str_return = lib.gillespie(param_input)
str_return2 = lib.gillespie(param_input)
str_return3 = lib.gillespie(param_input)
str_return4 = lib.gillespie(param_input)
str_return5 = lib.gillespie(param_input)
str_return6 = lib.gillespie(param_input)
str_return7 = lib.gillespie(param_input)
str_return8 = lib.gillespie(param_input)
print str_return.contents
print "*************"
print str_return2.contents
print "*************"
print str_return3.contents
print "*************"
print str_return4.contents
print "*************"
print str_return5.contents
print "*************"
print str_return6.contents
print "*************"
print str_return7.contents
print "*************"
print str_return8.contents
