import MooScripts as moo
import KimScripts as kim
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pylab
from scipy.stats.stats import pearsonr
from StringIO import StringIO

from ctypes import *
lib = cdll.LoadLibrary('./libfoo.so')
lib.gillespie.restype = c_char_p

# We pass in an i because over iterations the name of the dictionary changes from '0' to '1' ect.

# Simulate takes in a particle or set of particles and returns the y data that is simulated from them
# This data should be compared to the experimental x data to determine the correctness of any given particle
# Ref PNAS --> PRC2.1 "Generate a data set x** ~ f (x|theta)"
def simulate(param, i = 0):
    #print "BEGINNING SIMULATION OF:", param
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
    #print param_str
    #print theta_str
    param_num = theta_str.count(' ') + 1
    # Here we actually call the simulation script which skips through a python proxy to the C++ script one_trajectory
    #y = pd.read_table(StringIO(kim.ssa(param_str)), skiprows=0, nrows=5000, sep=" ", usecols = range(0,N_Species + 1))
    #print "Into hell"
    kim.ssa(param_str)
    #print "Out of hell"
    y = pd.read_csv("CtoPy.csv").ix[:,0:param_num]
    #print "Result is:", y

    # Establish a array for column names and ensure that the first one is Time because we will need to manually implement that column
    column_name = []
    column_name.append("Time")
    # For each species, make a new column title to mark how that specie changes over time.
    for x in range(0, N_Species):
        column_name.append('S'+str(x))
    y.columns = column_name
    y = y.set_index(['Time'])


    # y = y.iloc[0:-2] # Somehow the last time point becomes corrupted with NaN and other special characters.
    # make sure the bits that we want to be integer data type become that.
    for col in column_name[1:]:
        y[col] = y[col].apply(int)
    #print "ENDING SIMULATION, FORMATED RESULT IS:", y
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
    print "escape"
    return df
    # return {'Sinit':Sinit, 'theta':theta, 't_param':t_param}



    # param_inputs = pd.DataFrame()
    # for i in range(0, N_part):
    #     param_inputs.assign(i=pd.concat([Sinit, thetas[i], t_param]))
    # print param_inputs

# Here are warning variable which will be set to true if things go wrong but are corrected for.
WeightWarning = False

# math model is defined in c codes. 
# three species: init copy numbers = [10, 0, 100]
# 4 parameters: 
# example: species 1 synthesis rate, deg rate = 0.4, 2
# Second last param: time interval
# The last param: Number of data points.

# Simulation parameter initialization
# Here we manually generate an initial paramater particle style dataframe'
print "ONE TEST --------------------------------------------------------------"
#N_Species = 3
N_Species = 1
#N_param = 4
N_param = 2

# Define Each paramaters physical meaning when compared to the simulation program
#param_input1 = input_to_df("10 0 100 0.4 0.5 100 0.01 0.01 1750", N_Species, N_param)
param_input1 = input_to_df("10 120 0.8 0.025 120", N_Species, N_param)
#param_input1 = input_to_df("10 200 0.5 10 200", N_Species, N_param)

print "THREE TEST --------------------------------------------------------------"
#param_input3 = input_to_df('10 0 100 0.4 4 400 1 0.01 1750', N_Species, N_param)

param_input = param_input1

# abc algorithm optimization initialization
N_iter = 25
N_part = 20
Sims = 5
# Selects the number of slices in the time.
Slices = 10

# synthetic experimental data
print "TWO TEST --------------------------------------------------------------"
#param_input2 = input_to_df('10 0 100 0.4 1 200 0.01 0.01 1750', N_Species, N_param)
#param_input2 = input_to_df("10 100 1 0.025 120", N_Species, N_param)
#param_input2 = input_to_df("10 100 1 10 200", N_Species, N_param)

param_input2 = input_to_df("10 100 1 0.025 120", N_Species, N_param)

# tack the mean values of parameters among particles at given iteration points.
meanTrack = []
colorTrack = range(N_iter)
for colorDex in range(N_iter):
    colorTrack[colorDex] = 1 - (colorTrack[colorDex]/(N_iter*1.0))

print "simulate?"
x = simulate(param_input2)

#used for plotting the meanTrack, and for weighting in priorProb
#start = np.array([0.4, 0.5, 100, 0.1])
start = np.array([120, 0.8])
# Making the sigmas too small can result in null weight values which crashes the system.
sigmas = np.array([100, 700])
#truth = np.array([0.4, 1, 200, 0.01])
truth = np.array([100, 1])
print "BEGIN"

# real synthetic data generated (# = Sims)
xset = []
for xi in range(0, Sims):
    xset.append(simulate(param_input2))

# This generates the array of slice times for now
timeStep = xset[1].index.values[2] - xset[1].index.values[1]
timeMax = xset[1].index.values[-1]
SliceArray = [0]
for slicei in range(Slices - 1):
    SliceArray.append(SliceArray[slicei] + timeMax/Slices)
print xset[1].index.values
print SliceArray

# epsilon changes over iteration. 
# epsilon (iteration) is plotted.
# Now the PNAS process starts in earnest
# input: a threshold epsilon
# Ref PNAS --> PRC1
epsilon = 1500
epsilonTrack = [epsilon]
distanceArray = [-1]
failed = N_part
# Start the population indicator at 0 and let it run to the total number of allowed iterations
for t in range(0, N_iter):
    print "top of t=", t
    epsilon = kim.epsilon_next(epsilon, failed, N_part, distanceArray)
    epsilonTrack.append(epsilon)

    #print "t="
    #print t

    # on the first cycle, we will generate our initial particles from the starting particle
    # Ref PNAS --> PRC2.1 when t = 1
    if t == 0:
        #param_inputs = pd.Series([None] * N_part)
        param_inputs = kim.initial_particles(param_input, N_part)
        #print param_inputs
        #smack = plt.figure()
        #pylab.show()
        param_tilde = param_inputs.copy()
        paramSaved = param_inputs.copy()
        param = pd.DataFrame()
        # Set all weights equally? Seems a little out of place but could be a fence case
        w = pd.Series([1/float(N_part)]*N_part)
        wpast = w
        initParam = plt.figure()
        plt.plot(param_inputs.loc['theta'].loc['theta0'],param_inputs.loc['theta'].loc['theta1'],"ob")
        plt.plot(truth[0], truth[1], "or")
        plt.plot(start[0], start[1], "om")
        plt.xlabel("Paramater 0")
        plt.ylabel("Paramater 1")
        pylab.show()
        #print "INITIAL PARTICLES SELECTED"
        #print param_inputs
    # Otherwise we want to select some particles out of the previous set based on weights but with a random distribution
    # Then we want to use these selected particles to generate more particles using these selected according to the
    # perturbation kernel
    # Ref PNAS --> PRC2/1 when t > 0
    else:
        #print 't > 0'
        prevThet = param
        num_selected = int( 0.75*len(param.columns) )
        num_not_selected = len(param.columns) - num_selected
        #print "num_selected = ", num_selected
        # pick particles out of the bunch by weights in this case we don't pass actual particles in, just
        # the corresponding weights because the nature of the partcile does not affect the selection
        print "inputs"
        print w
        print num_selected
        param_input_index_selected = kim.select_particles(w, num_selected)
        #print "num of the selected particles = ", len(param_input_index_selected)
        print "param"
        print param
        temp = pd.DataFrame()
        # Now we build up temp, by appending each particle to it as they are selected by index
        print "?"
        print param_input_index_selected
        print zip(param_input_index_selected, range(0, len(param_input_index_selected)))
        for ind, k in zip(param_input_index_selected, range(0, len(param_input_index_selected))):
            #print "ind= ", ind, "k=", k
            temp[k] = param[ind]

        #print "Weights:"
        #print w
        print "PARTICLES SELECTED"
        #print temp

        # This block seperates the paramaters and initial values out of the particle dataframe from for temp
        print temp
        thetaTemp = temp.loc['theta']
        topSinit = temp.loc['Sinit']
        topTParam = temp.loc['t_param']
        thetaIndecies =  list(thetaTemp.index)

        # Now we generate new particles using the perturbation kernel
        # Number of particles not selected
        offset = (N_part - len(param_input_index_selected))
        #print temp
        #print "dtype?"
        #print "t = ", t
        
        # based on the 75% param sets, generate their statistics and generate 100% new params. 
        newThetas = pd.DataFrame(moo.PerturbationKernel(temp, N_part))
        newThetas.index = thetaIndecies
        #print newThetas

        for newThet in range(N_part):

            topSinit[newThet] = topSinit[1]
            topTParam[newThet] = topTParam[1]
            thetaTemp[newThet] = newThetas.iloc[:,newThet]

        perturbSave = temp.copy()
        print perturbSave
        print "Pause"
        #pause = plt.figure()
        #pylab.show()
        dict = {'Sinit': topSinit, 'theta': thetaTemp, 't_param': topTParam}
        temp = pd.concat(dict)
        #print "PERTURBATION COMPLETE, RESULT:"
        #print temp
        paramSaved = temp.copy()
        param_tilde = temp.copy()
        #param_inputs = perturbed according to a transition kernel K

    # Now that we have generated new particles we want to test each and every one of them and choose wheather or now they
    # be passed onto the next test. Here we set the particle indicator to 0 and begin to cycle through them as they appear
    # in temp or param?
    # Ref PNAS --> PRC2.1 Generate a data set and If p(S(x**)) both

    #print "param BEFORE IT ALL WAS UNDONE ----------------------------------------------"
    #print param
    #print "paramSaved"
    #print paramSaved
    #print "param_inputs?"
    #print param_inputs

    #print "NOW WE BEGIN THE CULLING BASED ON SIMULATION COMPARISON"
    i = 0
    # for each particle
    failed = 0
    distanceArray = []
    while (i < N_part):

        print "--------------------------"
        print "top of i=", i, "and t=", t

        # Run the simulation as detailed on the top of this python file and in the one_trajectory C++ script
        #print "PARTICLES WE HAVE NOW"
        #print param_tilde
        #print "param"
        #print param
        y = simulate(pd.DataFrame(param_tilde[i]), i)
        yset = []
        print "simulating..."
        for yi in range(0, Sims):
            yset.append(simulate(pd.DataFrame(param_tilde[i]), i))
        print "simulated"

        #print "yset"
        #print yset
        #print "y"
        #print y
        #print "DONE"
        #print yset[1]
        #print "SIMULATE OUT"
        #print y.head()
        # Normalize all the data points in x and y based on their max values. This is to remove any order of magnitude
        # issues between the expirimental and predicted data.
        x_norm = x.apply(moo.normalizeConditional)
        #print "xnorm"
        #print x_norm

        xset_norm = xset
        ##xset_norm = [bitx.apply(moo.normalizeConditional) for bitx in xset]
        #for xni in range(0, Sims):
            #xset_norm.append(xset[xni].apply(moo.normalizeConditional))

        y_norm = y.apply(moo.normalizeConditional)
        #print "ynorm"
        #print y_norm

        yset_norm = yset
        ##yset_norm = [bity.apply(moo.normalizeConditional) for bity in yset]

        #print "yset_norm"
        #print yset_norm
        #print "xset_norm"
        #print xset_norm
        print "t = ", t
        print "failed = ", failed

        #print "xset_norm[1]"
        #print xset_norm[1]
        #print "x_norm"
        #print x_norm
        num_elements = len(xset_norm[1].columns)*len(xset_norm[1].index)
        #num_elements = len(x_norm.columns)*len(x_norm.index)
        #print "num_elements"
        #print num_elements

        set_elements = len(xset_norm[1].columns)*len(xset_norm[1].index)
        #set_elements = len(x_norm.columns)*len(x_norm.index)
        #print "set_elements"
        #print set_elements

        #print "difference"
        #print (x_norm-y_norm)

        #exit(1)

        setDiffWhole = 0
        print "xset_norm[1] - yset_norm[1]"
        print xset_norm[1] - yset_norm[1]
        print "comparing..."
        for yyi in range(0, Sims):
            #for yyi in range(0, Sims):
            #tempDiff = [np.absolute(bitnx - yset_norm[yyi]) for bitnx in xset_norm]
            tempDiff = ((xset_norm[yyi] - yset_norm[yyi])**2).sum()
            print "tempDiff"
            print tempDiff
            #setDiffWhole = setDiffWhole + sum(bitcx.apply(lambda x:x**2).sum().sum() for bitcx in tempDiff)
            setDiffWhole = setDiffWhole + tempDiff
        print "compared."
                #setDiffWhole = setDiffWhole + (xset_norm[xxi]-yset_norm[xxi]).apply(lambda x:x**2).sum().sum()
        setDistance = setDiffWhole.iloc[0]**0.5
        # Set a new temp (Ovewrwriting the one generated above?) which represents the first step to the distance function
        # recording the errors
        distanceWhole = (x_norm-y_norm).apply(lambda x:x**2).sum().sum()
        #print "distanceWhole"
        #print distanceWhole

        #print "setDiffWhole"
        #print setDiffWhole

        print param_tilde[i]

        #distance = distanceWhole/float(num_elements)**0.5
        #setDistance = setDiffWhole/float(num_elements)**0.5
        print "within t=", t
        #print "distance for particle", i , "=", distance
        print "distance for set particle", i, "=", setDistance
        print "epsilon is,", epsilon
        print "particle"
        print param_tilde[i]

        print "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
        # This stuff is all for printing purposes.
        ysetLength = len(yset)
        print "ysetLength"
        print ysetLength


        #TimeLiney = y.index.values.tolist()
        TimeLiney = yset[1].index.values.tolist()
        Yzerosum = [0]*len(TimeLiney)
        fuckingThing = Yzerosum
        for yin in range(ysetLength):
            ybyte = yset[yin]
            Yzerobyte = ybyte['S0'].tolist()
            if (len(Yzerosum) > len(Yzerobyte)):
                fuckingThing = Yzerobyte
            for yin1 in range(len(fuckingThing)):
                Yzerosum[yin1] = Yzerosum[yin1] + Yzerobyte [yin1]
        Yzero = [float(yin2) / ysetLength for yin2 in Yzerosum]

        #Yone = y['S1']
        #Ytwo = y['S2']

        xsetLength = len(xset)
        print "xsetLength"
        print xsetLength
        #TimeLinex = x.index.values.tolist()
        TimeLinex = xset[1].index.values.tolist()
        Xzerosum = [0]*len(TimeLinex)
        for xin in range(xsetLength):
            xbyte = xset[xin]
            Xzerobyte = xbyte['S0'].tolist()
            for xin1 in range(len(Xzerosum)):
                Xzerosum[xin1] = Xzerosum[xin1] + Xzerobyte [xin1]
        Xzero = [float(xin2) / xsetLength for xin2 in Xzerosum]

        timespace = np.linspace(0, 3, num = 120)
        time0 = 10
        c = time0 - truth[0]/truth[1]

        amtspace = []
        for sindex in range(len(timespace)):
            bit = (c * np.exp(-truth[1]*timespace[sindex])) + truth[0]/truth[1]
            amtspace.append(bit)


        if t == 2:

            track = plt.figure()
            track1 = track.add_subplot(331)
            plt.plot(TimeLinex, Xzero, "r", TimeLiney, Yzero, "g", timespace, amtspace, "r--")
            plt.xlabel("Time")
            plt.ylabel("Concentration Species 1")
            plt.title("x/y trajectory")
            pylab.show()
        print "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"


        #if (t == 3):
            #exit(1)
        # If the distance it below tolerance its a good particle and we will save it to the param set.
        if np.absolute(setDistance) <= epsilon:

            print "PARTICLE", i, "PASSED"
            param[i] = param_tilde[i]
            if i==0:
                ys = {0:y}
            else:
                ys[i] = y
            distanceArray.append(np.absolute(setDistance))
            # Increase the value of the particle counter, if the particle fails the distance test then we are no closer
            # to a new generation of promising particles.
            i = i + 1
        # If not, we want to remove that particle and use the perturbation kernel to generate a new one from the same base
        # or re perturb the previous particle in a different random way in hopes of helping it pass
        else:

            print "PARTICLE", i, "FAILED"
            failed = failed + 1
            if failed > N_part*50:
                print "Particles failed exceeds particle number times 50, consider relaxing epsilon or speed of epsilon decrease."
                print param_tilde
                print param
                epsilon = epsilon * 1.025
                failed = 0
                #exit(1)
            #print i
            #print t
            #print "param_tilde"
            #print param_tilde
            if (t == 0):
                param_tilde[i] = kim.initial_particles(param_input, 1)
            else:
                # generate a new param_input
                # replace the below with a Kernel function
                #print "perturbSave"
                #print perturbSave
                newParams = moo.PerturbationKernel(perturbSave, 1)
                newParams = newParams[0]
                for j in range(len(newParams)):
                    param_tilde[i].loc['theta'][j] = newParams[j]
        print "-----------"
            #print "NEW PARTICLE"
            #print param_tilde[i]
            #print "when dist > epsilon, dummy, dummy"

    theta0i = param.loc['theta'].loc['theta0'].tolist()
    theta1i = param.loc['theta'].loc['theta1'].tolist()
    #theta2i = param.loc['theta'].loc['theta2'].tolist()
    #theta3i = param.loc['theta'].loc['theta3'].tolist()
    #mean = np.array([np.mean(theta0i), np.mean(theta1i), np.mean(theta2i), np.mean(theta3i)])
    mean = np.array([np.mean(theta0i), np.mean(theta1i)])
    meanTrack.append(mean)

    iterParam = plt.figure()
    plt.plot(param.loc['theta'].loc['theta0'],param.loc['theta'].loc['theta1'],"ob")
    plt.plot(truth[0], truth[1], "or")
    plt.plot(start[0], start[1], "om")
    plt.xlabel("Paramater 0")
    plt.ylabel("Paramater 1")
    # Comment out this if you want the script to stop showing a graph on each iteration
    pylab.show()

    #print "NOW WE ASSIGN NEW WEIGHTS"
    #print "CURRENT WEIGHTS"
    #print w
    #print "PASTWEIGHTS"
    #print wpast
    #calculate weight for all particles, if this is the first particle, then its weight is 1
    # I have backed this all down below the while loop... because both papers say we should update the thetas weights
    # all at once.
    wprep = w
    wnew = pd.Series([1/float(N_part)]*N_part)
    for k in range(N_part):
        if t == -1:
            #print "if"
            w = pd.Series([1/float(N_part)]*N_part)
            wpast = w
        # otherwise calculate the weight based on all of the existing weights, and all the past weights according to the
        # formula in the PNAS paper
        else:
            #print "top of k=", k
            #print "Weighting begins"
            #print "k"
            #print k
            #print "param"
            #print param
            #print "paramSaved"
            #print paramSaved
            #print "param_inputs?"
            #print param_inputs
            # Now we define weights based on prior weights, and all current shifts in weights and paramaters and the prior
            #print "N_part"
            #print N_part
            #print "wnew"
            #print wnew
            #print k
            wnew[k] = moo.weightt(N_part, w, wpast, param, paramSaved, k, t, start, sigmas)
            #print wnew
    print "particles failed culling"
    print failed

    # Normalize the wieghts
    print "Before norm"
    print wnew
    wnew = wnew/wnew.sum()
    print "wnew"
    print wnew
    if wnew.sum() == 0:
        print "WARNING: Weights near zero, sigma range must be widened or your prior set is too far off. Weights being reset to simple equal values"
        WeightWarning = True
        wnew = pd.Series([1/float(N_part)]*N_part)
        exit(1)

    #print "NEW WEIGHTS"
    #print wnew
    wpast = wprep
    w = wnew
# Having the first weight always be one massively corrupts the process resulting in the more iterations proceed
# the smaller the weights assigned to particles that are not the first one which results in the particles being picked
# almost always being 5 copies of the first particle.

# We still occasionally see corruption leaking into the program and crashing it? It happens somewhere in simulation
# but this means the tighter we make the system, the more runs it has to do, the greater the odds that it will crash.

# Any great number of iterations causes a certain particle to run away with weight and overwrite the other particles.
print param
print epsilon

print meanTrack

param.to_csv("ParticleTest.csv")
#exit(1)
print "0-1"
print pearsonr(param.loc['theta'].loc['theta0'].tolist(),param.loc['theta'].loc['theta1'].tolist())
#print "0-2"
#print pearsonr(param.loc['theta'].loc['theta0'].tolist(),param.loc['theta'].loc['theta2'].tolist())
#print "0-3"
#print pearsonr(param.loc['theta'].loc['theta0'].tolist(),param.loc['theta'].loc['theta3'].tolist())
#print "1-3"
#print pearsonr(param.loc['theta'].loc['theta1'].tolist(),param.loc['theta'].loc['theta3'].tolist())
#print "2-3"
#print pearsonr(param.loc['theta'].loc['theta2'].tolist(),param.loc['theta'].loc['theta3'].tolist())

theta0 = param.loc['theta'].loc['theta0'].tolist()
theta1 = param.loc['theta'].loc['theta1'].tolist()
#theta2 = param.loc['theta'].loc['theta2'].tolist()
#theta3 = param.loc['theta'].loc['theta3'].tolist()

means = np.array([np.mean(theta0), np.mean(theta1)]) #, np.mean(theta2), np.mean(theta3)])

ysetLength = len(yset)
print "ysetLength"
print ysetLength
#TimeLiney = y.index.values.tolist()
TimeLiney = yset[1].index.values.tolist()
Yzerosum = [0]*len(TimeLiney)
for yin in range(ysetLength):
    ybyte = yset[yin]
    Yzerobyte = ybyte['S0'].tolist()
    for yin1 in range(len(Yzerosum)):
        Yzerosum[yin1] = Yzerosum[yin1] + Yzerobyte [yin1]
Yzero = [float(yin2) / ysetLength for yin2 in Yzerosum]

#Yone = y['S1']
#Ytwo = y['S2']

xsetLength = len(xset)
print "xsetLength"
print xsetLength
#TimeLinex = x.index.values.tolist()
TimeLinex = xset[1].index.values.tolist()
Xzerosum = [0]*len(TimeLinex)
for xin in range(xsetLength):
    xbyte = xset[xin]
    Xzerobyte = xbyte['S0'].tolist()
    for xin1 in range(len(Xzerosum)):
        Xzerosum[xin1] = Xzerosum[xin1] + Xzerobyte [xin1]
Xzero = [float(xin2) / xsetLength for xin2 in Xzerosum]

print "YZERO"
print Yzero

print "Xzero"
print Xzero

#Xone = x['S1'].tolist()
#Xtwo = x['S2'].tolist()

timespace = np.linspace(0, 3, num = 120)
time0 = 10
c = time0 - truth[0]/truth[1]

amtspace = []
for i in range(len(timespace)):
    bit = (c * np.exp(-truth[1]*timespace[i])) + truth[0]/truth[1]
    amtspace.append(bit)

track = plt.figure()
track1 = track.add_subplot(331)
plt.plot(TimeLinex, Xzero, "r", TimeLiney, Yzero, "g", timespace, amtspace, "r--")
plt.xlabel("Time")
plt.ylabel("Concentration Species 1")
plt.title("x/y trajectory")

#track2 = track.add_subplot(332)
#plt.plot(TimeLine, Xone, "g")
#plt.xlabel("Time")
#plt.ylabel("Concentration Species 2")
#plt.title("One_Trajectory")

#track3 = track.add_subplot(333)
#plt.plot(TimeLine, Xtwo, "r")
#plt.xlabel("Time")
#plt.ylabel("Concentration Species 3")
#plt.title("One_Trajectory")

oneToTwo = track.add_subplot(339)
plt.plot(param.loc['theta'].loc['theta0'],param.loc['theta'].loc['theta1'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[0], meanBit[1], "o", color = str(colorTrack[meani]))
plt.plot(truth[0], truth[1], "or")
plt.plot(means[0], means[1], "og")
plt.plot(start[0], start[1], "om")
plt.xlabel("Paramater 0")
plt.ylabel("Paramater 1")
#plt.title("P0 vs P1")

pylab.show()

stats = plt.figure()
plt.plot(range(N_iter + 1), epsilonTrack, "-r")
plt.xlabel("iterations")
plt.ylabel("epsilon value")
pylab.show()

exit(1)

oneToThree = track.add_subplot(335)
plt.plot(param.loc['theta'].loc['theta0'],param.loc['theta'].loc['theta2'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[0], meanBit[2], "o", color = str(colorTrack[meani]))
plt.plot(truth[0], truth[2], "or")
plt.plot(means[0], means[2], "og")
plt.plot(start[0], start[2], "om")
plt.xlabel("Paramater 0")
plt.ylabel("Paramater 2")
#plt.title("P0 vs P2")

oneToFour = track.add_subplot(336)
plt.plot(param.loc['theta'].loc['theta0'],param.loc['theta'].loc['theta3'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[0], meanBit[3], "o", color = str(colorTrack[meani]))
plt.plot(truth[0], truth[3], "or")
plt.plot(means[0], means[3], "og")
plt.plot(start[0], start[3], "om")

plt.xlabel("Paramater 0")
plt.ylabel("Paramater 3")
#plt.title("P0 vs P3")

twoToThree = track.add_subplot(337)
plt.plot(param.loc['theta'].loc['theta1'],param.loc['theta'].loc['theta2'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[1], meanBit[2], "o", color = str(colorTrack[meani]))
plt.plot(truth[1], truth[2], "or")
plt.plot(means[1], means[2], "og")
plt.plot(start[1], start[2], "om")
plt.xlabel("Paramater 1")
plt.ylabel("Paramater 2")
#plt.title("P1 vs P2")

twoToFour = track.add_subplot(338)
plt.plot(param.loc['theta'].loc['theta1'],param.loc['theta'].loc['theta3'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[1], meanBit[3], "o", color = str(colorTrack[meani]))
plt.plot(truth[1], truth[3], "or")
plt.plot(means[1], means[3], "og")
plt.plot(start[1], start[3], "om")
plt.xlabel("Paramater 1")
plt.ylabel("Paramater 3")
#plt.title("P1 vs P3")


threeToFour = track.add_subplot(339)
plt.plot(param.loc['theta'].loc['theta2'],param.loc['theta'].loc['theta3'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[2], meanBit[3], "o", color = str(colorTrack[meani]))
plt.plot(truth[2], truth[3], "or")
plt.plot(means[2], means[3], "og")
plt.plot(start[2], start[3], "om")
plt.xlabel("Paramater 2")
plt.ylabel("Paramater 3")
#plt.title("P2 vs P3")

print "Blue = Generated Particles"
print "Green = Generated Average"
print "Grayscale = Generated averages over time white -> black"
print "Magenta = Initial Particles"
print "Red = True Values"

pylab.show()

stats = plt.figure()
plt.plot(range(N_iter + 1), epsilonTrack, "-r")
plt.xlabel("iterations")
plt.ylabel("epsilon value")
pylab.show()

print param

if WeightWarning:
    print "WARNING: Weights near zero, sigma range must be widened or your prior set is too far off. Weights being reset to simple equal values"