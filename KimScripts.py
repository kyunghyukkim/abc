#!/Anaconda/python
__author__ = 'Kyung Kim'
import numpy

###############
# generate time series data from stochastic simulation algorithm
###############
from ctypes import *
lib = cdll.LoadLibrary('./libfoo.so')
lib.gillespie.restype = c_char_p

def ssa(param_input):
	# math model is defined in c codes. 
	# three species: init copy numbers = [10, 0, 100]
	# 4 parameters: 
	# example: species 1 synthesis rate, deg rate = 0.4, 2
	# Third last param: time interval
	# Second last param: Number of data points.
	# the last param: not necessary at this point. 
	
	str_return = lib.gillespie(param_input)
	return str_return
