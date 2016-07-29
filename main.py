import MooScripts as moo
import KimScripts as kim
import pandas as pd
from StringIO import StringIO

from ctypes import *
lib = cdll.LoadLibrary('./libfoo.so')
lib.gillespie.restype = c_char_p

def simulate(param, Num_Species):
    y = pd.read_table(StringIO(kim.ssa(param)), skiprows=0, nrows=5000, sep=" ", usecols = range(0,Num_Species+1))

    column_name = []
    column_name.append("Time")
    for x in range(0, Num_Species):
        column_name.append('S'+str(x))
    y.columns = column_name
    y = y.set_index(['Time'])
    return y

def initialize_particle_dist(param_input, N_Species, N_param, N_part):
    # sample particles with given initial values of theta (obtained from param_input).
    # kim.initial_particles(initial_theta, sample_size)
    param_input_list = map(float, param_input.split(' '))
    S_delt = param_input_list[0:N_Species + 1]
    S_delt = map(int, param_input_list[0:N_Species]) + [float(param_input_list[N_Species])]
    print "Initial Species copy numbers and time step used in the simulation algorithm\n", S_delt
    theta = param_input_list[N_Species + 1:N_Species + N_param + 1]
    print "Initial theta value = ", theta
    runtime = param_input_list[N_Species + N_param + 1:]
    # print runtime

    thetas = kim.initial_particles(theta, N_part)
    param_inputs = [None] * N_part
    for i in range(0, N_part):
        lst = S_delt + thetas[i].tolist() + runtime
        param_inputs[i] = ' '.join(map(str, lst))
    print "First set of particles"
    print param_inputs
    return param_inputs

# math model is defined in c codes. 
# three species: init copy numbers = [10, 0, 100]
# 4 parameters: 
# example: species 1 synthesis rate, deg rate = 0.4, 2
# Third last param: time interval
# Second last param: Number of data points.
# the last param: not necessary at this point. 

N_Species = 3
N_param = 4
param_input = '10 0 100 0.4 2 100 0.1 0.02 10 240000'
N_iter = 1
N_part = 3

# synthetic experimental data
param_input2 = '10 0 100 0.4 1 200 0.01 0.2 10 240000'
x = simulate(param_input2, N_Species)
print x

# input: a threshold epsilon
epsilon = 0.1

for t in range(1, N_iter+1):
    epsilon = kim.epsilon_next(epsilon)

    param_inputs = [None] * N_part
    if t == 1:
        param_inputs = initialize_particle_dist(param_input, N_Species, N_param, N_part)
    else:
        print 't > 1'#param_inputs = from previous population with weights
        #param_inputs = perturbed according to a transition kernel K
    i = 0
    # while (i < N_part):
    #     #print param_inputs[i]
    #     y = simulate(param_inputs[i], N_Species)
    #     print y.head()
    #     print (x.unstack()-y.unstack()).apply(lambda x: x**2).sum()
    #
    #     #if moo.euclidd(pd.Series([1, 2, 3]), pd.Series([4, 3, 1]))
    #
    #     i = i + 1
    #     #print type(y)
    #     #print  y.rstrip('\0')



#    theta = kim.kernel(theta)

# #	determine thhe parameters of the perturbation kernel K(.|.)
# ############
# ##  This is justa rough draft.
# ##  Keagan, can you make the rest of the code?
# ############
# #	i=1
# #	repeat
# #		if t=1:
# 			sample theta_tilde from pi(theta)
# 		else:
# 			sample theta from the rpevious population with weights
# 			sample theta_tilde from K(.|theta) and such that pi(theta_tilde) > 0
# 		sample y from f(.|theta_tilde)
# 		if moo.euclidd(y,x) <= epsilon:
# 			theta[i][t] = theat_tilde
# 			y[i][t] = y
# # stochastic simulation
# # something like this...
# tseries =  kim.ssa(param_input)
#
