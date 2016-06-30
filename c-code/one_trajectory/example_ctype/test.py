from ctypes import * 
lib = cdll.LoadLibrary('./libfoo.so')

#f.bar(c_double(1.2))
lib.f.restype = c_float
print lib.f(c_float(1.3))
