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
    y = y.iloc[0:-2] # Somehow the last time point becomes corrupted with NaN and other special characters.
    for col in column_name[1:]:
        y[col] = y[col].apply(int)
    return y

def initialize_particle_dist(param_input, N_Species, N_param, N_part):
    # sample particles with given initial values of theta (obtained from param_input).
    # kim.initial_particles(initial_theta, sample_size)
    param_input_list = map(float, param_input.split(' '))
    Sinit = map(int, param_input_list[0:N_Species])
    print "Initial Species copy numbers used in the simulation algorithm\n", Sinit
    theta = param_input_list[N_Species:N_Species + N_param]
    print "Initial theta value = ", theta
    delt = [param_input_list[N_Species + N_param]]
    print "Simulation time step = ", delt
    runtime = map(int, param_input_list[N_Species + N_param + 1:])
    # print runtime

    thetas = kim.initial_particles(theta, N_part)
    param_inputs = pd.Series([None] * N_part)
    for i in range(0, N_part):
        lst = Sinit + thetas[i].tolist() + delt + runtime
        param_inputs[i] = ' '.join(map(str, lst))
    #print "First set of particles"
    #print param_inputs
    return param_inputs, thetas



# math model is defined in c codes. 
# three species: init copy numbers = [10, 0, 100]
# 4 parameters: 
# example: species 1 synthesis rate, deg rate = 0.4, 2
# Second last param: time interval
# The last param: Number of data points.


N_Species = 3
N_param = 4
param_input = '10 0 100 0.4 2 100 0.1 0.02 10'
N_iter = 2
N_part = 5

# synthetic experimental data
param_input2 = '10 0 100 0.4 1 200 0.01 0.02 10'
x = simulate(param_input2, N_Species)

# input: a threshold epsilon
epsilon = 1.0

for t in range(0, N_iter):
    epsilon = kim.epsilon_next(epsilon)


    if t == 0:
        param_inputs = pd.Series([None] * N_part)
        param_inputs, thetas_tilde = initialize_particle_dist(param_input, N_Species, N_param, N_part)
        thetas = pd.DataFrame(0, columns = thetas_tilde.columns, index = thetas_tilde.index)
        w = pd.Series([0.0]*N_part)
        print "t=0"
        # print param_inputs
    else:
        print 't > 0'
        num_selected = int( 0.75*len(thetas.columns) )
        num_not_selected = len(thetas.columns) - num_selected
        print "num_selected = ", num_selected
        param_input_index_selected = kim.select_particles(thetas, w, num_selected)
        print "num of the selected particles = ", len(param_input_index_selected)

        temp = pd.Series()
        print param_input_index_selected
        print param_inputs
        for ind, k in zip(param_input_index_selected, range(0, len(param_input_index_selected))):
            #print "ind= ", ind, "k=", k
            #print param_inputs[ind]
            temp.set_value(k, param_inputs[ind])
        for k in range(len(param_input_index_selected), N_part):
            exit(1)
            #temp.set_value(k, get_param_from_K())
        #print param_inputs
        #print param_inputs[{0,0,1}]
        exit(1)
        #param_inputs = perturbed according to a transition kernel K


    i = 0
    while (i < N_part):
        print "value of i = ", i, "input = ", param_inputs[i]
        y = simulate(param_inputs[i], N_Species)
        print y.head()

        # distance computation
        x_norm = x.apply(lambda x: x*len(x)/x.sum())
        y_norm = y.apply(lambda y: y*len(y)/y.sum())
        num_elements = len(x_norm.columns)*len(x_norm.index)
        temp = (x_norm-y_norm).apply(lambda x:x**2).sum().sum()
        distance = temp/float(num_elements)**0.5
        print "distance", distance
        if distance <= epsilon:
            thetas[i] = thetas_tilde[i]
            if i==0:
                ys = {0:y}
            else:
                ys[i] = y
            i = i + 1
        else:
            # generate a new param_input
            # replace the below with a Kernel function
            # [param_inputs[i]], [thetas_tilde[i]] = initialize_particle_dist(param_input, N_Species, N_param, 1)
            print "when dist > epsilon, dummy, dummy"
            # print param_inputs
        #calculate weight for all particles
        if i == 0:
            w[i] = 1
        else:
            w[i] = 1+2
        w = w/w.sum()