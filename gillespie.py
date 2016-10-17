from ctypes import *

# This line causes the error (Keagan)
lib = cdll.LoadLibrary('./libfoo.so')

##f.bar(c_double(1.2))
#lib.f.restype = c_float
#print lib.f(c_float(1.3))

lib.gillespie.restype = c_char_p
# math model is defined in c codes. 
# three species: init copy numbers = [10, 0, 100]
# 4 parameters: 
# example: species 1 synthesis rate, deg rate = 0.4, 2
# Third last param: time interval
# Second last param: Number of data points.
# the last param: not necessary at this point. 

param_input = '10 0  100  0.4    2   100   0.1  0.02   10   240000'
str_return = lib.gillespie(param_input)
str_return2 = lib.gillespie(param_input)
print str_return
print "*************"
print str_return2
