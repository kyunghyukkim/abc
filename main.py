import MooScripts as moo
import KimScripts as kim
import pandas as pd
from StringIO import StringIO

from ctypes import *
lib = cdll.LoadLibrary('./libfoo.so')
lib.gillespie.restype = c_char_p

def simulate(param):
    print "SIMULATE ----------------------------------------------------------"
    N_Species = len(param.loc['Sinit'])
    N_param = len(param.loc['theta'])

    Sinit_str = ' '.join(map(str, param.loc['Sinit'][0].tolist()))
    theta_str = ' '.join(map(str, param.loc['theta'][0].tolist()))
    t_param_str = ' '.join(map(str, param.loc['t_param'][0].tolist()))

    param_str = ' '.join([Sinit_str, theta_str, t_param_str])
    print param_str
    print "SIM INTER ---------------------------------------------------------"
    y = pd.read_table(StringIO(kim.ssa(param_str)), skiprows=0, nrows=5000, sep=" ", usecols = range(0,N_Species + 1))
    print "MARK --------------------------------------------------------------"
    print y
    column_name = []
    column_name.append("Time")
    print "SIM IN2ER ---------------------------------------------------------"
    for x in range(0, N_Species):
        column_name.append('S'+str(x))
    print "SIM IN3ER ---------------------------------------------------------"
    print "Bluh"
    print y
    y.columns = column_name
    print "SIMCOL"
    y = y.set_index(['Time'])
    print y
    print "SIMDEX"
    y = y.iloc[0:-2] # Somehow the last time point becomes corrupted with NaN and other special characters.
    print y
    print "SIM IN4ER ---------------------------------------------------------"
    for col in column_name[1:]:
        print y[col]
        y[col] = y[col].apply(int)
    print "SIM IN5ER ---------------------------------------------------------"
    return y



def input_to_df(input_str, N_Species, N_param):
    input_list = map(float, input_str.split(' '))
    theta_index_name = []
    for x in range(0, N_param):
        theta_index_name.append('theta' + str(x))
    Sinit_index_name = []
    for x in range(0, N_Species):
        Sinit_index_name.append('S' + str(x))

    Sinit = pd.DataFrame(map(int, input_list[0:N_Species]), index = Sinit_index_name, dtype='object')
    print "Initial Species copy numbers used in the simulation algorithm\n", Sinit

    theta = pd.DataFrame(input_list[N_Species:N_Species + N_param], index = theta_index_name, dtype='object')
    print "Initial theta values\n", theta

    t_param = pd.DataFrame(map(float, [input_list[N_Species + N_param]]), index = ['delt'], dtype='object')
    t_param.loc['N_time']=map(int, [input_list[N_Species + N_param+1]])
    print "Simulation time parameters\n", t_param
    dict = {'Sinit': Sinit, 'theta': theta, 't_param': t_param}
    df = pd.concat(dict)
    return df
    # return {'Sinit':Sinit, 'theta':theta, 't_param':t_param}



    # param_inputs = pd.DataFrame()
    # for i in range(0, N_part):
    #     param_inputs.assign(i=pd.concat([Sinit, thetas[i], t_param]))
    # print param_inputs



# math model is defined in c codes. 
# three species: init copy numbers = [10, 0, 100]
# 4 parameters: 
# example: species 1 synthesis rate, deg rate = 0.4, 2
# Second last param: time interval
# The last param: Number of data points.

# Simulation parameter initialization
print "ONE TEST --------------------------------------------------------------"
N_Species = 3
N_param = 4
param_input = input_to_df("10 0 100 0.4 2 100 0.1 0.02 10", N_Species, N_param)


# abc algorithm optimization initialization
N_iter = 2
N_part = 5

# synthetic experimental data
print "TWO TEST --------------------------------------------------------------"
param_input2 = input_to_df('10 0 100 0.4 1 200 0.01 0.02 10', N_Species, N_param)
print "SIM TEST --------------------------------------------------------------"
x = simulate(param_input2)


# input: a threshold epsilon
epsilon = 1.0

print "FOR TEST --------------------------------------------------------------"
print range(0, N_iter)
print "Format Check"
print param_input2
for t in range(0, N_iter):
    epsilon = kim.epsilon_next(epsilon)
    print "t="
    print t

    if t == 0:
        #param_inputs = pd.Series([None] * N_part)
        param_inputs = kim.initial_particles(param_input, N_part)
        param_tilde = param_inputs
        param = pd.DataFrame()
        w = pd.Series([0.0]*N_part)
        print "t=0"
        # print param_inputs
    else:
        print 't > 0'
        num_selected = int( 0.75*len(param.columns) )
        num_not_selected = len(param.columns) - num_selected
        print "num_selected = ", num_selected
        param_input_index_selected = kim.select_particles(w, num_selected)
        print "num of the selected particles = ", len(param_input_index_selected)

        temp = pd.DataFrame()
        print param_input_index_selected
        print param
        for ind, k in zip(param_input_index_selected, range(0, len(param_input_index_selected))):
            #print "ind= ", ind, "k=", k
            #print param_inputs[ind]
            temp[k] = param[ind]
        #print temp
        temp.join(get_params_from_K(N_part - lon(param_input_index_selected)))
        #exit(1)
        #param_inputs = perturbed according to a transition kernel K

    print "WHILE TEST --------------------------------------------------------------"
    i = 0
    print "N+part"
    print N_part
    exit(1)
    while (i < N_part):
        print "i="
        print i
        print "TILDE TES~T"
        print pd.DataFrame(param_tilde[0])
        #print "RENAME TEST"
        #print pd.DataFrame(param_tilde[1].rename(0))
        print "SIMULATE IN"
        y = simulate(pd.DataFrame(param_tilde[i]))
        print "SIMULATE OUT"
        print y.head()
        # distance computation
        x_norm = x.apply(lambda x: x*len(x)/x.sum())
        y_norm = y.apply(lambda y: y*len(y)/y.sum())
        num_elements = len(x_norm.columns)*len(x_norm.index)
        temp = (x_norm-y_norm).apply(lambda x:x**2).sum().sum()
        distance = temp/float(num_elements)**0.5
        print "distance", distance
        if distance <= epsilon:
            param[i] = param_tilde[i]
            if i==0:
                ys = {0:y}
            else:
                ys[i] = y
            i = i + 1
        else:
            # generate a new param_input
            # replace the below with a Kernel function
            param_tilde[i] = kim.initial_particles(param_input, 1)[0]
            #print "when dist > epsilon, dummy, dummy"
            # print param_inputs
        #calculate weight for all particles
        if i == 0:
            w[i] = 1
        else:
            w[i] = 1+2
        w = w/w.sum()
