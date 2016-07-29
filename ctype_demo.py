# from ctypes import *
# lib = cdll.LoadLibrary('./libfoo.so')
#
# ##f.bar(c_double(1.2))
# #lib.f.restype = c_float
# #print lib.f(c_float(1.3))
#
# lib.gillespie.restype = c_char_p
# # math model is defined in c codes.
# # three species: init copy numbers = [10, 0, 100]
# # 4 parameters:
# # example: species 1 synthesis rate, deg rate = 0.4, 2
# # Third last param: time interval
# # Second last param: Number of data points.
# # the last param: not necessary at this point.
#
# param_input = '10 0  100  0.4    2   100   0.1  0.02   10   240000'
# str_return = lib.gillespie(param_input)
# #str_return2 = lib.gillespie(param_input)
# print "STR_OUT = ", str_return
#

# #!/usr/bin/python2.7
# -*- coding: utf-8 -*-
 
"""Basic ctypes demo."""
 
import ctypes
import ctypes.util

mylib = ctypes.cdll.LoadLibrary('libmything.dylib')
 
libc = ctypes.cdll.LoadLibrary(ctypes.util.find_library('c'))

param_input = '10 0  100  0.4    2   100   0.1  0.02   10   240000'

# Configure return type of say_hi to "pointer to char"
# mylib.say_hi.restype = ctypes.POINTER(ctypes.c_char)
# hi_msg_p = mylib.say_hi('OOgi')
#
mylib.gillespie.restype = ctypes.POINTER(ctypes.c_char)
sim_output = mylib.gillespie(param_input)
# Get the string from the char-pointer
sim_output = ctypes.string_at(sim_output)
print "second method"
print sim_output
print "second method after"
# libc.free(sim_output)

sim_output2 = mylib.gillespie(param_input)
# Get the string from the char-pointer
sim_output2 = ctypes.string_at(sim_output2)
print "second call"
print sim_output2
print "second call after"
# cleanup
