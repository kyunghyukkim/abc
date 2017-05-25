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
    param_num = theta_str.count(' ') + 1
    specie_num = Sinit_str.count(' ') + 2 # plus 2 for the inbetween space and the time column
    # Here we actually call the simulation script which skips through a python proxy to the C++ script one_trajectory
    #If you want to return to the direct C++ to python method use this line. y = pd.read_table(StringIO(kim.ssa(param_str)), skiprows=0, nrows=5000, sep=" ", usecols = range(0,N_Species + 1))
    kim.ssa(param_str)
    y = pd.read_csv("CtoPy.csv").ix[:,0:specie_num]

    # Establish a array for column names and ensure that the first one is Time because we will need to manually implement that column
    column_name = []
    column_name.append("Time")
    # For each species, make a new column title to mark how that specie changes over time.
    for x in range(0, N_Species):
        column_name.append('S'+str(x))
    y.columns = column_name
    y = y.set_index(['Time'])


    # make sure the bits that we want to be integer data type become that.
    for col in column_name[1:]:
        y[col] = y[col].apply(int)
    return y


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

# Here are warning variable which will be set to true if things go wrong but are corrected for.
WeightWarning = False

# math model is defined in c codes. 
# three species: init copy numbers = [10, 0, 100]
# 4 parameters: 
# example: species 1 synthesis rate, deg rate = 0.4, 2
# Second last param: time interval
# The last param: Number of data points.

#########################################################################################################################################################################################################
                                                                                                                                                                                                        #
# Simulation parameter initialization                                                                                                                                                                   #
# Here we manually generate an initial paramater particle style dataframe'                                                                                                                              #
N_Species = 6                                                                                                                                                                                           #
N_param = 4                                                                                                                                                                                             #
                                                                                                                                                                                                        #
#For Graphing, set to be equal to your initial guess for each parameter                                                                                                                                 #
start = np.array([70, 30, 60, 30])                                                                                                                                                                      #
#For perturbation and error calculations, making these too small can cause a crash which will return a specific error to let you know, an order of magnitude greater than you think is probably safe.   #
sigmas = np.array([10000, 100, 100, 100])                                                                                                                                                               #
#For Graphing if you don't have experimental data set these to be equal to the true values of your parameters                                                                                           #
truth = np.array([200, 10, 500, 4])                                                                                                                                                                     #
                                                                                                                                                                                                        #
                                                                                                                                                                                                        #
# smc algorithm optimization initialization                                                                                                                                                             #
# Number of iterations for the whole process, recommended 10                                                                                                                                            #
N_iter = 3                                                                                                                                                                                              #
# Number of particles generated and maintained, recommended 50                                                                                                                                          #
N_part = 50                                                                                                                                                                                             #
# Number of simulations run to compare data sets, recommended 10                                                                                                                                        #
Sims = 10                                                                                                                                                                                               #
# Initial epsilon to be reduced from, changes drastically depends on error calculation method and circuit                                                                                               #
epsilon = 2500                                                                                                                                                                                          #
# Selects the number of slices in the time. Currently inactivated for ineffectiveness.                                                                                                                  #
Slices = 4                                                                                                                                                                                              #
                                                                                                                                                                                                        #
#########################################################################################################################################################################################################

# Define Each paramaters physical meaning when compared to the simulation program, this defines the guess you'd like to make
#format([Initial Conc], [Priors], [steps, number of steps])
print "TESTTTT"
print "0 0 0 0 0 0 " + str(start[0]) + " " + str(start[1]) + " " + str(start[2]) + " " + str(start[3])+ " " + "0.00025 2000"
print "0 0 0 0 0 0 70 30 60 30 0.00025 2000"
param_input1 = input_to_df("0 0 0 0 0 0 " + str(start[0]) + " " + str(start[1]) + " " + str(start[2]) + " " + str(start[3])+ " " + "0.00025 2000", N_Species, N_param)
param_input = param_input1

# Defines a bound on set graphs because of the random nature of the gillespie algorithms.
complen = param_input1.loc['t_param'][0].loc['N_time'] - 1

# synthetic experimental data
print "TWO TESTSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSs"
param_input2 = input_to_df("0 0 0 0 0 0 " + str(truth[0]) + " " + str(truth[1]) + " " + str(truth[2]) + " " + str(truth[3])+ " " + "0.00025 2000", N_Species, N_param)
# track the mean values of parameters among particles at given iteration points.
meanTrack = []
colorTrack = range(N_iter)
for colorDex in range(N_iter):
    colorTrack[colorDex] = 1 - (colorTrack[colorDex]/(N_iter*1.0))
x = simulate(param_input2)

# Generates the synthetic experimental data
xset = []
for xi in range(0, Sims):
    xset.append(simulate(param_input2))


# This generates the array of slice times for now
timeStep = xset[1].index.values[2] - xset[1].index.values[1]
timeMax = xset[1].index.values[-1]
SliceArray = [0]
for slicei in range(Slices):
    SliceArray.append(SliceArray[slicei] + timeMax/Slices)

# real synthetic data generated (# = Sims)
# Generate all of the test data, place the slices into sets and place those sets into a list.a
xset = []
xSlicesSet = []
endpointsSet = []
for xi in range(0, Sims):
    xSlices = []
    endpoints = []
    simBit = simulate(param_input2)
    xset.append(simBit)
    for si in range(Slices):
        xSlices.append(simBit[SliceArray[si]: SliceArray[si + 1]])
    xSlicesSet.append(xSlices)


# epsilon changes over iteration.
# epsilon (iteration) is plotted.
# Now the PNAS process starts in earnest
# input: a threshold epsilon
# Ref PNAS --> PRC1

epsilonTrack = []
distanceArray = [-1]
failed = N_part
# Start the population indicator at 0 and let it run to the total number of allowed iterations
for t in range(0, N_iter):
    print "top of iteration t=", t
    epsilon = kim.epsilon_next(epsilon, failed, N_part, distanceArray)
    epsilonTrack.append(epsilon)

    # on the first cycle, we will generate our initial particles from the starting particle
    # Ref PNAS --> PRC2.1 when t = 1
    if t == 0:

        param_inputs = kim.initial_particles(param_input, N_part)

        param_tilde = param_inputs.copy()
        paramSaved = param_inputs.copy()
        param = pd.DataFrame()
        # Set all weights equally? Seems a little out of place but could be a fence case
        w = pd.Series([1/float(N_part)]*N_part)
        wpast = w

        # Plot initial parameter distribution
        initParam = plt.figure()
        init01 = initParam.add_subplot(321)
        plt.plot(param_inputs.loc['theta'].loc['theta0'],param_inputs.loc['theta'].loc['theta1'],"ob")
        plt.plot(truth[0], truth[1], "or")
        plt.plot(start[0], start[1], "om")
        plt.xlabel("alpha")
        plt.ylabel("k")

        init01 = initParam.add_subplot(322)
        plt.plot(param_inputs.loc['theta'].loc['theta0'],param_inputs.loc['theta'].loc['theta2'],"ob")
        plt.plot(truth[0], truth[2], "or")
        plt.plot(start[0], start[2], "om")
        plt.xlabel("alpha")
        plt.ylabel("beta")

        init01 = initParam.add_subplot(323)
        plt.plot(param_inputs.loc['theta'].loc['theta0'],param_inputs.loc['theta'].loc['theta3'],"ob")
        plt.plot(truth[0], truth[3], "or")
        plt.plot(start[0], start[3], "om")
        plt.xlabel("alpha")
        plt.ylabel("n")

        init01 = initParam.add_subplot(324)
        plt.plot(param_inputs.loc['theta'].loc['theta1'],param_inputs.loc['theta'].loc['theta2'],"ob")
        plt.plot(truth[1], truth[2], "or")
        plt.plot(start[1], start[2], "om")
        plt.xlabel("k")
        plt.ylabel("beta")

        init01 = initParam.add_subplot(325)
        plt.plot(param_inputs.loc['theta'].loc['theta1'],param_inputs.loc['theta'].loc['theta3'],"ob")
        plt.plot(truth[1], truth[3], "or")
        plt.plot(start[1], start[3], "om")
        plt.xlabel("k")
        plt.ylabel("n")

        init01 = initParam.add_subplot(326)
        plt.plot(param_inputs.loc['theta'].loc['theta2'],param_inputs.loc['theta'].loc['theta3'],"ob")
        plt.plot(truth[2], truth[3], "or")
        plt.plot(start[2], start[3], "om")
        plt.xlabel("beta")
        plt.ylabel("n")

        pylab.show()


    # Otherwise we want to select some particles out of the previous set based on weights but with a random distribution
    # Then we want to use these selected particles to generate more particles using these selected according to the
    # perturbation kernel
    # Ref PNAS --> PRC2/1 when t > 0
    else:
        prevThet = param.copy()
        num_selected = int( 0.75*len(param.columns) )
        num_not_selected = len(param.columns) - num_selected

        # pick particles out of the bunch by weights in this case we don't pass actual particles in, just
        # the corresponding weights because the nature of the partcile does not affect the selection
        param_input_index_selected = kim.select_particles(w, num_selected)

        temp = pd.DataFrame()
        # Now we build up temp, by appending each particle to it as they are selected by index
        for ind, k in zip(param_input_index_selected, range(0, len(param_input_index_selected))):
            temp[k] = param[ind]

        print "PARTICLES SELECTED"

        # This block seperates the paramaters and initial values out of the particle dataframe from for temp
        print temp
        thetaTemp = temp.loc['theta']
        topSinit = temp.loc['Sinit']
        topTParam = temp.loc['t_param']
        thetaIndecies = list(thetaTemp.index)

        # Now we generate new particles using the perturbation kernel
        # Number of particles not selected
        offset = (N_part - len(param_input_index_selected))

        # based on the 75% param sets, generate their statistics and generate 100% new params. 
        newThetas = pd.DataFrame(moo.PerturbationKernel(temp, N_part))
        newThetas.index = thetaIndecies

        for newThet in range(N_part):

            topSinit[newThet] = topSinit[1]
            topTParam[newThet] = topTParam[1]
            thetaTemp[newThet] = newThetas.iloc[:,newThet]

        perturbSave = temp.copy()
        dict = {'Sinit': topSinit, 'theta': thetaTemp, 't_param': topTParam}
        temp = pd.concat(dict)
        paramSaved = temp.copy()
        param_tilde = temp.copy()

    # Now that we have generated new particles we want to test each and every one of them and choose wheather or now they
    # be passed onto the next test. Here we set the particle indicator to 0 and begin to cycle through them as they appear
    # in temp or param?
    #print "NOW WE BEGIN THE CULLING BASED ON SIMULATION COMPARISON"
    i = 0
    # for each particle
    failed = 0
    distanceArray = []
    while (i < N_part):

        print "--------------------------"
        print "top of i=", i, "and t=", t

        # Run the simulation as detailed on the top of this python file and in the one_trajectory C++ script

        print SliceArray
        print xSlicesSet[0][0].loc[:,'S3'].iloc[0]

        SliceByte = param_input2.loc['t_param'][0][1]/Slices
        yset = []
        ySlicesSet = []
        for yi in range(0, Sims):
            ySlices = []
            simBit = simulate(pd.DataFrame(param_tilde[i]), i)
            yset.append(simBit)
            # Still must be hand edited for each system.
            for si in range(Slices):
                param_inputnew = param_tilde[i].copy()
                param_inputnew.loc['t_param'][1] = SliceByte
                param_inputnew.loc['Sinit'][0] = xSlicesSet[yi][si].loc[:,'S0'].iloc[0]
                param_inputnew.loc['Sinit'][1] = xSlicesSet[yi][si].loc[:,'S1'].iloc[0]
                param_inputnew.loc['Sinit'][2] = xSlicesSet[yi][si].loc[:,'S2'].iloc[0]
                param_inputnew.loc['Sinit'][3] = xSlicesSet[yi][si].loc[:,'S3'].iloc[0]
                param_inputnew.loc['Sinit'][4] = xSlicesSet[yi][si].loc[:,'S4'].iloc[0]
                param_inputnew.loc['Sinit'][5] = xSlicesSet[yi][si].loc[:,'S5'].iloc[0]
                ySlices.append(simulate(pd.DataFrame(param_inputnew), i))
            ySlicesSet.append(ySlices)


        y = simulate(pd.DataFrame(param_tilde[i]), i)
        yset = []
        for yi in range(0, Sims):
            yset.append(simulate(pd.DataFrame(param_tilde[i]), i))

        # Normalize all the data points in x and y based on their max values. This is to remove any order of magnitude
        # issues between the expirimental and predicted data.
        x_norm = x#.apply(moo.normalizeConditional)
        #print "xnorm"
        #print x_norm

        xset_norm = xset
        ##xset_norm = [bitx.apply(moo.normalizeConditional) for bitx in xset]
        #for xni in range(0, Sims):
            #xset_norm.append(xset[xni].apply(moo.normalizeConditional))

        y_norm = y#.apply(moo.normalizeConditional)

        yset_norm = yset
        ##yset_norm = [bity.apply(moo.normalizeConditional) for bity in yset]

        print "t = ", t
        print "failed = ", failed

        num_elements = len(xset_norm[1].columns)*len(xset_norm[1].index)

        set_elements = len(xset_norm[1].columns)*len(xset_norm[1].index)

        sliceDiffWhole = 0
        for syi in range(Sims):
            for slyi in range(Slices - 1):
                tempDiff = ((xSlicesSet[syi][slyi].loc[:,'S3'].values[:SliceByte-1] - ySlicesSet[syi][slyi].loc[:,'S3'].values[:SliceByte-1])**2).sum()
                print "temp diff for", syi, "and", slyi, "=", tempDiff
                sliceDiffWhole = sliceDiffWhole + tempDiff
        sliceDistance = sliceDiffWhole**0.5


        setDiffWhole = 0
        for yyi in range(0, Sims):
            tempDiff = ((xset_norm[yyi] - yset_norm[yyi])**2).sum()
            setDiffWhole = setDiffWhole + tempDiff
        setDistance = setDiffWhole.iloc[0]**0.5
        # Set a new temp (Ovewrwriting the one generated above?) which represents the first step to the distance function
        # recording the errors
        distanceWhole = (x_norm-y_norm).apply(lambda x:x**2).sum().sum()

        #print "param_tilde"
        #print param_tilde
        #plt.figure()
        #pylab.show()

        #distance = distanceWhole/float(num_elements)**0.5
        #setDistance = setDiffWhole/float(num_elements)**0.5
        print "within t=", t
        #print "distance for particle", i , "=", distance
        print "distance for set particle", i, "=", setDistance#sliceDistance
        print "epsilon is,", epsilon
        print "particle"
        print param_tilde[i]
        '''
        print "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
        ysetLength = len(ySlicesSet)
        print "ysetLength"
        print ysetLength

        YStichedSet = []
        XStichedSet = []
        lengthOfOdd = len(ySlicesSet[0][0]) * Slices * Slices
        for pti in range(Sims):
            YPreStitch = ySlicesSet[pti]
            XPreStitch = xSlicesSet[pti]
            YStitchStart = YPreStitch[0]
            XStitchStart = XPreStitch[0]


            for sti in range(Slices)[1:]:
                XStitchStart = XStitchStart.append(XPreStitch[sti])
                YStitchStart = YStitchStart.append(YPreStitch[sti])
            YStichedSet.append(YStitchStart)
            XStichedSet.append(XStitchStart)
            prevOdd = lengthOfOdd
            lengthOfOdd = min(len(YStitchStart), len(XStitchStart), prevOdd)



        #TimeLiney = y.index.values.tolist()
        print lengthOfOdd

        TimeLiney = YStitchStart.index.values.tolist()
        TimeLinex = XStitchStart.index.values.tolist()

        Yzerosum = [0]*lengthOfOdd

        for yin in range(ysetLength):
            ybyte = YStichedSet[yin]
            Yzerobyte = ybyte['S3'].tolist()[:complen]
            for yin1 in range(len(Yzerosum)):
                Yzerosum[yin1] = Yzerosum[yin1] + Yzerobyte[yin1]
        Yzero = [float(yin2) / ysetLength for yin2 in Yzerosum]

        #Yone = y['S1']
        #Ytwo = y['S2']

        xsetLength = len(xSlicesSet)
        print "xsetLength"
        print xsetLength
        #TimeLinex = x.index.values.tolist()

        Xzerosum = [0]*complen
        for xin in range(xsetLength):
            xbyte = XStichedSet[xin]
            Xzerobyte = xbyte['S3'].tolist()[:complen]
            for xin1 in range(len(Yzerosum)):
                Xzerosum[xin1] = Xzerosum[xin1] + Xzerobyte [xin1]
        Xzero = [float(xin2) / xsetLength for xin2 in Xzerosum]

        print "YZERO"
        print Yzero

        print "Xzero"
        print Xzero

        #Xone = x['S1'].tolist()
        #Xtwo = x['S2'].tolist()

        track = plt.figure()
        track1 = track.add_subplot(331)
        plt.plot(TimeLinex[:lengthOfOdd], Xzero[:lengthOfOdd], "r", TimeLinex[:lengthOfOdd], Yzero[:lengthOfOdd], "g")
        plt.xlabel("Time")
        plt.ylabel("Concentration Species 1")
        plt.title("x/y trajectory")
        print "setDistance"
        print setDistance
        print "Particle:"
        print param_tilde[i]
        pylab.show()
        print "|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
        '''
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
    theta2i = param.loc['theta'].loc['theta2'].tolist()
    theta3i = param.loc['theta'].loc['theta3'].tolist()
    mean = np.array([np.mean(theta0i), np.mean(theta1i), np.mean(theta2i), np.mean(theta3i)])
    #mean = np.array([np.mean(theta0i), np.mean(theta1i)])
    meanTrack.append(mean)

    iterParam = plt.figure()
    iter01 = iterParam.add_subplot(321)
    plt.plot(param.loc['theta'].loc['theta0'],param.loc['theta'].loc['theta1'],"ob")
    plt.plot(truth[0], truth[1], "or")
    plt.plot(start[0], start[1], "om")
    plt.xlabel("Paramater 0")
    plt.ylabel("Paramater 1")

    iter01 = iterParam.add_subplot(322)
    plt.plot(param.loc['theta'].loc['theta0'],param.loc['theta'].loc['theta2'],"ob")
    plt.plot(truth[0], truth[2], "or")
    plt.plot(start[0], start[2], "om")
    plt.xlabel("Paramater 0")
    plt.ylabel("Paramater 2")

    iter01 = iterParam.add_subplot(323)
    plt.plot(param.loc['theta'].loc['theta0'],param.loc['theta'].loc['theta3'],"ob")
    plt.plot(truth[0], truth[3], "or")
    plt.plot(start[0], start[3], "om")
    plt.xlabel("Paramater 0")
    plt.ylabel("Paramater 3")

    iter01 = iterParam.add_subplot(324)
    plt.plot(param.loc['theta'].loc['theta1'],param.loc['theta'].loc['theta2'],"ob")
    plt.plot(truth[1], truth[2], "or")
    plt.plot(start[1], start[2], "om")
    plt.xlabel("Paramater 1")
    plt.ylabel("Paramater 2")

    iter01 = iterParam.add_subplot(325)
    plt.plot(param.loc['theta'].loc['theta1'],param.loc['theta'].loc['theta3'],"ob")
    plt.plot(truth[1], truth[3], "or")
    plt.plot(start[1], start[3], "om")
    plt.xlabel("Paramater 1")
    plt.ylabel("Paramater 3")

    iter01 = iterParam.add_subplot(326)
    plt.plot(param.loc['theta'].loc['theta2'],param.loc['theta'].loc['theta3'],"ob")
    plt.plot(truth[2], truth[3], "or")
    plt.plot(start[2], start[3], "om")
    plt.xlabel("Paramater 2")
    plt.ylabel("Paramater 3")

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
        if t == N_iter - 1:
            #print "if"
            w = pd.Series([1/float(N_part)]*N_part)
            wpast = w
        # otherwise calculate the weight based on all of the existing weights, and all the past weights according to the
        # formula in the PNAS paper
        else:

            wnew[k] = moo.weightt(N_part, w, wpast, param, paramSaved, k, t, start, sigmas)
            #print wnew
    print "particles failed culling"
    print failed

    # Normalize the wieghts
    wnew = wnew/wnew.sum()
    print "New Weights"
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
theta2 = param.loc['theta'].loc['theta2'].tolist()
theta3 = param.loc['theta'].loc['theta3'].tolist()

means = np.array([np.mean(theta0), np.mean(theta1), np.mean(theta2), np.mean(theta3)])


ysetLength = len(ySlicesSet)
print "ysetLength"
print ysetLength

YStichedSet = []
XStichedSet = []
lengthOfOdd = len(ySlicesSet[0][0]) * Slices * Slices
for pti in range(Sims):
    YPreStitch = ySlicesSet[pti]
    XPreStitch = xSlicesSet[pti]
    YStitchStart = YPreStitch[0]
    XStitchStart = XPreStitch[0]


    for sti in range(Slices)[1:]:
        XStitchStart = XStitchStart.append(XPreStitch[sti])
        YStitchStart = YStitchStart.append(YPreStitch[sti])
    YStichedSet.append(YStitchStart)
    XStichedSet.append(XStitchStart)
    prevOdd = lengthOfOdd
    lengthOfOdd = min(len(YStitchStart), len(XStitchStart), prevOdd)



#TimeLiney = y.index.values.tolist()
print lengthOfOdd

TimeLiney = YStitchStart.index.values.tolist()
TimeLinex = XStitchStart.index.values.tolist()

Yzerosum = [0]*lengthOfOdd

for yin in range(ysetLength):
    ybyte = YStichedSet[yin]
    Yzerobyte = ybyte['S3'].tolist()[:complen]
    for yin1 in range(len(Yzerosum)):
        Yzerosum[yin1] = Yzerosum[yin1] + Yzerobyte[yin1]
Yzero = [float(yin2) / ysetLength for yin2 in Yzerosum]

xsetLength = len(xSlicesSet)
print "xsetLength"
print xsetLength
#TimeLinex = x.index.values.tolist()

Xzerosum = [0]*complen
for xin in range(xsetLength):
    xbyte = XStichedSet[xin]
    Xzerobyte = xbyte['S3'].tolist()[:complen]
    for xin1 in range(len(Yzerosum)):
        Xzerosum[xin1] = Xzerosum[xin1] + Xzerobyte [xin1]
Xzero = [float(xin2) / xsetLength for xin2 in Xzerosum]

print "YZERO"
print Yzero

print "Xzero"
print Xzero

track = plt.figure()
track1 = track.add_subplot(331)
plt.plot(TimeLinex[:lengthOfOdd], Xzero[:lengthOfOdd], "r", TimeLinex[:lengthOfOdd], Yzero[:lengthOfOdd], "g")
plt.xlabel("Time")
plt.ylabel("Concentration Species 1")
plt.title("x/y trajectory")
pylab.show()

print "SET-----------------------------------------------------"

ysetLength = len(yset)
print "ysetLength"
print ysetLength
#TimeLiney = y.index.values.tolist()


TimeLiney = yset[1].index.values.tolist()
TimeLinex = xset[1].index.values.tolist()


Yzerosum = [0]*complen
for yin in range(ysetLength):
    ybyte = yset[yin]
    Yzerobyte = ybyte['S3'].tolist()
    for yin1 in range(len(Yzerosum)):
        Yzerosum[yin1] = Yzerosum[yin1] + Yzerobyte [yin1]
Yzero = [float(yin2) / ysetLength for yin2 in Yzerosum]

#Yone = y['S1']
#Ytwo = y['S2']

xsetLength = len(xset)
print "xsetLength"
print xsetLength
#TimeLinex = x.index.values.tolist()

Xzerosum = [0]*complen
for xin in range(xsetLength):
    xbyte = xset[xin]
    Xzerobyte = xbyte['S3'].tolist()
    for xin1 in range(len(Yzerosum)):
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

#Defunct for repressilator
amtspace = []
for i in range(len(timespace)):
    bit = (c * np.exp(-truth[1]*timespace[i])) + truth[0]/truth[1]
    amtspace.append(bit)

track = plt.figure()
track1 = track.add_subplot(333)
plt.plot(TimeLinex[:complen], Xzero[:complen], "r", TimeLiney[:complen], Yzero[:complen], "g") #, timespace, amtspace, "r--")
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
print "meanTrack"
print meanTrack
print meanTrack[0]
print meanTrack[0][2]

zeroToOne = track.add_subplot(334)
plt.plot(param.loc['theta'].loc['theta0'],param.loc['theta'].loc['theta1'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[0], meanBit[1], "o", color = str(colorTrack[meani]))
plt.plot(truth[0], truth[1], "or")
plt.plot(means[0], means[1], "og")
plt.plot(start[0], start[1], "om")
plt.xlabel("alpha")
plt.ylabel("k")

zeroToTwo = track.add_subplot(335)
plt.plot(param.loc['theta'].loc['theta0'],param.loc['theta'].loc['theta2'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[0], meanBit[2], "o", color = str(colorTrack[meani]))
plt.plot(truth[0], truth[2], "or")
plt.plot(means[0], means[2], "og")
plt.plot(start[0], start[2], "om")
plt.xlabel("alpha")
plt.ylabel("beta")

zeroToThree = track.add_subplot(336)
plt.plot(param.loc['theta'].loc['theta0'],param.loc['theta'].loc['theta3'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[0], meanBit[3], "o", color = str(colorTrack[meani]))
plt.plot(truth[0], truth[3], "or")
plt.plot(means[0], means[3], "og")
plt.plot(start[0], start[3], "om")
plt.xlabel("alpha")
plt.ylabel("n")

oneToTwo = track.add_subplot(337)
plt.plot(param.loc['theta'].loc['theta1'],param.loc['theta'].loc['theta2'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[1], meanBit[2], "o", color = str(colorTrack[meani]))
plt.plot(truth[1], truth[2], "or")
plt.plot(means[1], means[2], "og")
plt.plot(start[1], start[2], "om")
plt.xlabel("k")
plt.ylabel("beta")

oneToThree = track.add_subplot(338)
plt.plot(param.loc['theta'].loc['theta1'],param.loc['theta'].loc['theta3'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[1], meanBit[3], "o", color = str(colorTrack[meani]))
plt.plot(truth[1], truth[3], "or")
plt.plot(means[1], means[3], "og")
plt.plot(start[1], start[3], "om")
plt.xlabel("k")
plt.ylabel("n")

twoToThree = track.add_subplot(339)
plt.plot(param.loc['theta'].loc['theta2'],param.loc['theta'].loc['theta3'],"ob")
for meani in range(N_iter):
    meanBit = meanTrack[meani]
    plt.plot(meanBit[2], meanBit[3], "o", color = str(colorTrack[meani]))
plt.plot(truth[2], truth[3], "or")
plt.plot(means[2], means[3], "og")
plt.plot(start[2], start[3], "om")
plt.xlabel("beta")
plt.ylabel("n")

pylab.show()

stats = plt.figure()
plt.plot(range(N_iter), epsilonTrack, "-r")
plt.xlabel("iterations")
plt.ylabel("epsilon value")
pylab.show()

print param

if WeightWarning:
    print "WARNING: Weights near zero, sigma range must be widened or your prior set is too far off. Weights being reset to simple equal values"