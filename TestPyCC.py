__author__ = 'Keagan Moo'

from ctypes import *

lib = cdll.LoadLibrary('./Testfoo.so')
chartype = c_char_p

lib.TestVoid.restype = chartype
lib.TestVoid.argtypes = [chartype]

parinput = "Wat?"

lib.TestVoid(parinput)

lib.TestStr.restype = chartype
lib.TestStr.argtypes = [chartype]

parinput = "Wat?"

result = lib.TestStr(parinput)
print(result)

