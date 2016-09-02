__author__ = 'Keagan Moo'

from ctypes import *

lib = cdll.LoadLibrary('./Testfoo.so')
chartype = c_char

lib.TestPrint()

