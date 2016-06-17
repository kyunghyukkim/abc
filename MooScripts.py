#!/Anaconda/python
__author__ = 'Keagan Moo'
import numpy

# euclidd (Euclidean Distance Function) will take two lists of points
# and calculate the distance between the two using the euclidian distance
# function. Which is described by sqrt((x1 - x2)^2+(y1 - y2)^2 ...)
def euclidd(list1, list2):
    iter = len(list1)
    sum = 0
    for i in range(iter):
        gap = (list1[i] - list2[i]) ** 2
        sum = sum + gap
    return sum ** (0.5)

# weightt (Weight Calculator) takes a tolerance counter value,
# a total number of particles, a previous weight vector,
# a theta distribution, the previous theta distribution,
# a perturbation function, and a proposal distribuiton to
# create a new weighting vector for particles.
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

# Kpert (K perturbation function) stands in for the current perturbation kernel
# it is defined online, but could be changed if needed. It takes the total number
# of particles, a standard deviation and the current and previous parameter distributions
# (theta).
def Kpert(N, sigj, pThet, pThetp):
    product = 1
    for j in range(N):
        bit = numpy.exp(-((pThet[j] - pThetp[j])**2)/((2*sigj)**2))/numpy.sqrt(2*3.14159*sigj)
        product = product * bit
    return product

# sigmaj produces the standard deviation for the Kpert function from the various statistical
# properties of the system.
def sigmaj(var):
    return var