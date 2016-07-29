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
    # print str_return
    return str_return


from math import *
import itertools
from matplotlib.pylab import *
import scipy.stats
import pandas as pd


"""
initialize the particles of theta with uniformly distributed random number between passed_on_value/10 and
passed_on_value*10.
If the passed_on_value is 0, then we assign random numbers between 10^-6 and 10^-5.

initial_particles( dimension of theta, number of particles)

For example, initial_particles([0,10,100,4], 10).
This returns pandas DataFrame structure as shown below.
                  0            1            2            3           4  \
theta0  6.50465e-06  1.30545e-06  6.00645e-06  1.51863e-06  5.8941e-06
theta1      1.11194      77.6596      47.1341      60.5072     15.1167
theta2      40.5187      587.044      456.305      601.678     859.439
theta3      33.3591      25.8456      20.9166       31.972     13.0722

                  5            6            7            8            9
theta0  6.65912e-06  4.58307e-06  4.95579e-06  4.43802e-06  3.56429e-06
theta1      67.2929      21.1594      34.6431        60.97      46.8012
theta2      87.2928      380.012      539.675      360.102      967.428
theta3      30.9787      28.9047      36.9889      38.9364      30.5391
"""
def initial_particles(theta, n):
    index_name = []
    for x in range(0, len(theta)):
        index_name.append('theta'+str(x))
    df = pd.DataFrame(columns=range(0,n), index=index_name)
    for i,index in enumerate(index_name):
    #    print i, index
        for col in range(0,n):
            if theta[i] > 0:
                df.set_value(index, col, numpy.random.uniform(theta[i]/10.,theta[i]*10,1)[0])
            else:
                df.set_value(index, col, numpy.random.uniform(10**-6,10**-5,1)[0])
    #print df
    return df

#initial_particles([0,10,100,4], 10)

def epsilon_next(epsilon):
    return epsilon

def kernel(theta):

    return new_theta



def sampler(n):
    """
    This function samples from x and returns a vetor shorter that len(vector size)
    """
    x = numpy.random.uniform(0, 1, n)  # Sample
    y = numpy.random.uniform(0, 1, n) * 1.5  # Envelope
    fx = scipy.stats.norm.pdf(x)  # Target. x has a beta distribution

    #print fx
    #plt.hist(fx)
    #show()
    #s = itertools.compress(x, y<fx)
    #print s # return only those values that satisfy the condition
    #return s


def efficiency(vector, n):
    """
    tests the efficiency of the sampling procedure.
    returns acceptance probability
    """
    l = len(vector)
    prob = l / float(n)
    diff = n - l
    n2 = int(diff / prob)  # n required to obtain the remaining samples needed
    vec2 = sampler(n2)
    s = concatenate((vector, vec2))
    # generates histogram
    nb, bins, patches = hist(s, bins=50, normed=0)

    xlabel('s')
    ylabel('frequency')
    title('Histogram of s: n=%s' % n)
    show()
    return s


n = 100000
sample = sampler(n)
print sample
# efficiency(sample, n)
