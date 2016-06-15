#!/Anaconda/python
__author__ = 'Keagan Moo'
import numpy


def euclidd(list1, list2):
    iter = len(list1)
    sum = 0
    for i in range(iter):
        gap = (list1[i] - list2[i]) ** 2
        sum = sum + gap
    return sum ** (0.5)

def weightt(t, N, wgtp, pThet, pThetp, pTurb, piTurb):
    if t == 0:
        return 1
    elif t > 0:
        sum = 0
        for i in range(N):
            bit = wgtp * pTurb(pThetp, pThet)
            sum = sum + bit
        return piTurb(pThet)/sum
    else:
        return -666

def Kpert(N, sigj, pThet, pThetp):
    product = 1
    for j in range(N):
        bit = numpy.exp(-((pThet[j] - pThetp[j])**2)/((2*sigj)**2))/numpy.sqrt(2*3.14159*sigj)
        product = product * bit
    return product


def sigmaj(var):
    return var