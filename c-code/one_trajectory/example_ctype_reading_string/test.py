from ctypes import * 
lib = cdll.LoadLibrary('./libfoo.so')

##f.bar(c_double(1.2))
#lib.f.restype = c_float
#print lib.f(c_float(1.3))

lib.f.restype = c_char_p
str = '10, 0, 100, 0,4, 2, 100, 0.1, 0.02, 10000, 240000'
str_return = lib.f(str)
print [float(x) for x in  str_return.split(',')]
#lst = [x  for x in lib.f('1, 2, 3')]
#print lst
