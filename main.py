import MooScripts as moo
import KimScripts as kim
import pandas as pd
from StringIO import StringIO

from ctypes import *
lib = cdll.LoadLibrary('./libfoo.so')
lib.gillespie.restype = c_char_p

# We pass in an i because over iterations the name of the dictionary changes from '0' to '1' ect.

# Simulate takes in a particle or set of particles and returns the y data that is simulated from them
# This data should be compared to the experimental x data to determine the correctness of any given particle
# Ref PNAS --> PRC2.1 "Generate a data set x** ~ f (x|theta)"
def simulate(param, i = 0):
    print "BEGINNING SIMULATION OF:", param
    # N_Species is the number of total species being considered. Sinit is the intial value of each specie
    N_Species = len(param.loc['Sinit'])
    # N_param is the number of parameters we are tracking each of which is evaluated in the theta part of the particle
    N_param = len(param.loc['theta'])

    # This block of code takes each portion of the Particles DataFrame and puts it in a coresponding variable
    Sinit_str = ' '.join(map(str, param.loc['Sinit'][i].tolist()))
    theta_str = ' '.join(map(str, param.loc['theta'][i].tolist()))
    t_param_str = ' '.join(map(str, param.loc['t_param'][i].tolist()))

    # param_str is now created as a single string representation of a single particle
    param_str = ' '.join([Sinit_str, theta_str, t_param_str])
    print param_str

    # Here we actually call the simulation script which skips through a python proxy to the C++ script one_trajectory
    y = pd.read_table(StringIO(kim.ssa(param_str)), skiprows=0, nrows=5000, sep=" ", usecols = range(0,N_Species + 1))
    print "Result is:", y

    # Establish a array for column names and ensure that the first one is Time because we will need to manually implement that column
    column_name = []
    column_name.append("Time")
    # For each species, make a new column title to mark how that specie changes over time.
    for x in range(0, N_Species):
        column_name.append('S'+str(x))
    y.columns = column_name
    y = y.set_index(['Time'])


    y = y.iloc[0:-2] # Somehow the last time point becomes corrupted with NaN and other special characters.
    # make sure the bits that we want to be integer data type become that.
    for col in column_name[1:]:
        y[col] = y[col].apply(int)
    print "ENDING SIMULATION, FORMATED RESULT IS:", y
    return y


# Revisit if convergance doesn't work. Change on the distibution of param inputs selected initially from prior input

# input to df is a formatting function that takes in the input string and the number of species and parameters for reference
# to turn the string form particle into a fully formated particle dataframe using dataframes within dataframes.
def input_to_df(input_str, N_Species, N_param):
    # input_str is the string representation of a particle, a list of initial conditions and parameters.
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
# Here we manually generate an initial paramater particle style dataframe
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

print "THREE TEST --------------------------------------------------------------"
param_input3 = input_to_df('10 0 100 0.4 4 400 1 0.02 10', N_Species, N_param)

x = simulate(param_input2)

# Now the PNAS process starts in earnest
# input: a threshold epsilon
# Ref PNAS --> PRC1
epsilon = 1

# Start the population indicator at 0 and let it run to the total number of allowed iterations
for t in range(0, N_iter):

    # later this function will update the tolerance to a lower value, for now it returns the same again
    epsilon = kim.epsilon_next(epsilon)
    print "t="
    print t

    # on the first cycle, we will generate our initial particles from the starting particle
    # Ref PNAS --> PRC2.1 when t = 1
    if t == 0:
        #param_inputs = pd.Series([None] * N_part)
        param_inputs = kim.initial_particles(param_input3, N_part)
        param_tilde = param_inputs
        paramSaved = param_inputs
        param = pd.DataFrame()
        # Set all weights equally? Seems a little out of place but could be a fence case
        w = pd.Series([1/float(N_part)]*N_part)
        print "t=0"
        print "INITIAL PARTICLES SELECTED"
        print param_inputs
        # print param_inputs
    # Otherwise we want to select some particles out of the previous set based on weights but with a random distribution
    # Then we want to use these selected particles to generate more particles using these selected according to the
    # perturbation kernel
    # Ref PNAS --> PRC2/1 when t > 0
    else:
        print 't > 0'
        prevThet = param
        num_selected = int( 0.75*len(param.columns) )
        num_not_selected = len(param.columns) - num_selected
        print "num_selected = ", num_selected
        # pick particles out of the bunch by weights in this case we don't pass actual particles in, just
        # the corresponding weights because the nature of the partcile does not affect the selection
        param_input_index_selected = kim.select_particles(w, num_selected)
        print "num of the selected particles = ", len(param_input_index_selected)

        temp = pd.DataFrame()
        # Now we build up temp, by appending each particle to it as they are selected by index
        for ind, k in zip(param_input_index_selected, range(0, len(param_input_index_selected))):
            #print "ind= ", ind, "k=", k
            #print param_inputs[ind]
            temp[k] = param[ind]

        print "Weights:"
        print w
        print "PARTICLES SELECTED"
        print temp

        # This block seperates the paramaters and initial values out of the particle dataframe from for temp
        thetaTemp = temp.loc['theta']
        topSinit = temp.loc['Sinit']
        topTParam = temp.loc['t_param']
        thetaIndecies =  list(thetaTemp.index)

        # Now we generate new particles using the perturbation kernel (DO WE WANT TO USE PARAM OR TEMP?)
        newThetas = pd.DataFrame(moo.PerturbationKernel(temp, N_part - len(param_input_index_selected)))
        newThetas.index = thetaIndecies

        for newThet in range(N_part - len(param_input_index_selected)):

            topSinit[newThet + N_part - 2] = topSinit[newThet]
            topTParam[newThet + N_part - 2] = topTParam[newThet]
            thetaTemp[newThet + N_part - 2] = newThetas.iloc[:,newThet]


        print temp
        dict = {'Sinit': topSinit, 'theta': thetaTemp, 't_param': topTParam}
        temp = pd.concat(dict)
        print "PERTURBATION COMPLETE, RESULT:"
        print temp
        paramSaved = temp
        #param_inputs = perturbed according to a transition kernel K

    # Now that we have generated new particles we want to test each and every one of them and choose wheather or now they
    # be passed onto the next test. Here we set the particle indicator to 0 and begin to cycle through them as they appear
    # in temp or param?
    # Ref PNAS --> PRC2.1 Generate a data set and If p(S(x**)) both
    print "NOW WE BEGIN THE CULLING BASED ON SIMULATION COMPARISON"
    i = 0
    wpast = w
    # for each particle
    while (i < N_part):
        print "i="
        print i
        # Run the simulation as detailed on the top of this python file and in the one_trajectory C++ script
        print "PARTICLES WE HAVE NOW"
        print param_tilde
        y = simulate(pd.DataFrame(param_tilde[i]), i)
        print "SIMULATE OUT"
        print y.head()
        # Normalize all the data points in x and y based on their max values. This is to remove any order of magnitude
        # issues between the expirimental and predicted data.
        x_norm = x.apply(lambda x: x*len(x)/x.sum())
        y_norm = y.apply(lambda y: y*len(y)/y.sum())
        num_elements = len(x_norm.columns)*len(x_norm.index)

        # Set a new temp (Ovewrwriting the one generated above?) which represents the first step to the distance function
        # recording the errors
        temp = (x_norm-y_norm).apply(lambda x:x**2).sum().sum()
        distance = temp/float(num_elements)**0.5
        print "distance for particle", i , "=", distance
        print "epsilon is,", epsilon

        # If the distance it below tolerance its a good particle and we will save it to the param set.
        if distance <= epsilon:
            param[i] = param_tilde[i]
            if i==0:
                ys = {0:y}
            else:
                ys[i] = y
            # Increase the value of the particle counter, if the particle fails the distance test then we are no closer
            # to a new generation of promising particles.
            i = i + 1
        # If not, we want to remove that particle and use the perturbation kernel to generate a new one from the same base
        # or re perturb the previous particle in a different random way in hopes of helping it pass
        else:

            print "PARTICLE", i, "FAILED"
            print i
            print t

            if (t == 0):
                param_tilde[i] = kim.initial_particles(param_input, 1)
            else:
                # generate a new param_input
                # replace the below with a Kernel function
                newParams = moo.PerturbationKernel(paramSaved, 1)
                newParams = newParams[0]
                for j in range(len(newParams)):
                    param_tilde[i].loc['theta'][j] = newParams[j]
            print "NEW PERTURBED PARTICLE"
            print param_tilde[i]
            #print "when dist > epsilon, dummy, dummy"
        # After generating a new particle... or passing a sucessful one, we generate a new weight for it.

    print "NOW WE ASSIGN NEW WEIGHTS"
    print "CURRENT WEIGHTS"
    print w
    #calculate weight for all particles, if this is the first particle, then its weight is 1
    # I have backed this all down below the while loop... because both papers say we should update the thetas weights
    # all at once.
    wprep = w
    for k in range(len(param.loc['theta'])):
        # Why is the first weight always 1?
        if k == -1:
            print "if"
            w[k] = 1
        # otherwise calculate the weight based on all of the existing weights, and all the past weights according to the
        # formula in the PNAS paper
        else:
            print "Weighting begins"
            print "k"
            print k

            # Now we define weights based on prior weights, and all current shifts in weights and paramaters and the prior
            print "N_part"
            print N_part
            print "wpast"
            w[k] = moo.weightt(N_part, w, wpast, param, paramSaved, k)

    # Normalize the wieghts
    w = w/w.sum()
    print "NEW WEIGHTS"
    print w
    wpast = wprep
# Having the first weight always be one massively corrupts the process resulting in the more iterations proceed
# the smaller the weights assigned to particles that are not the first one which results in the particles being picked
# almost always being 5 copies of the first particle.

# We still occasionally see corruption leaking into the program and crashing it? It happens somewhere in simulation
# but this means the tighter we make the system, the more runs it has to do, the greater the odds that it will crash.

# Any great number of iterations causes a certain particle to run away with weight and overwrite the other particles.
print param
param.to_csv("ParticleTest.csv")