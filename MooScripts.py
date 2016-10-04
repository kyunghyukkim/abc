#!/Anaconda/python
__author__ = 'Keagan Moo'
import numpy
import pandas

# euclidd (Euclidean Distance Function) will take two lists of points
# Kim: to speed up the computation by using vetors, lists need to be pandas.Series.
# and calculate the distance between the two using the euclidian distance
# function. Which is described by sqrt((x1 - x2)^2+(y1 - y2)^2 ...)
def euclidd(list1, list2):
    iter = len(list1)
    # sum = 0
    # for i in range(iter):
    #     gap = (list1[i] - list2[i]) ** 2
    #     sum = sum + gap
    # return sum ** (0.5)
    diff = list1 - list2
    return (diff.apply(lambda x: x**2).sum())**0.5


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

# Change to multi-variate normal
def Kpert(t, N, sigj, pThet, pThetp):
    product = 1
    for j in range(N):
        bit = numpy.exp(-((pThet[j] - pThetp[j])**2)/((2*sigj)**2))/numpy.sqrt(2*3.14159*sigj)
        product = product * bit
    return product

def PerturbationKernel(particles, N):

    print "Analyzing input type..."
    print type(particles)

    print "ViewingParticles..."
    print particles

    print "Extraction..."
    ThetaSet = particles.loc['theta']

    print "Viewing ThetaSet..."
    print ThetaSet
    print "Analyzing Dimensions..."
    rowCount = ThetaSet.shape[0]
    colCount = ThetaSet.shape[1]
    print "Viewing rowCount..."
    print rowCount
    print colCount
    means = []
    index = 0
    while index < rowCount:
        subjectParam = ThetaSet.iloc[[index]]
        print "Viewing subjectParam..."
        print subjectParam

        print "Appending Means..."
        currentMean = numpy.mean(subjectParam.as_matrix())
        print currentMean
        means = numpy.append(means, currentMean)
        print "Viewing Means..."
        print means

        print "Increment"
        index += 1
    print "Final Index..."
    print index
    print "Expected Index..."
    print rowCount
        #CovarPrep = pandas.DataFrame([CovarPrep,subjectParam])
        #print "Building Covars..."
        #CovarPrep[label] = subjectParam
    print "Building Covars..."
    Covar = numpy.cov(ThetaSet)
    print "Viewing Covars..."
    print Covar

    print "Generating Results..."
    windex = 0
    results = pandas.DataFrame()
    while windex < N:

        result = numpy.random.multivariate_normal(means, Covar)
        results[windex] = result
        windex += 1

    print "Viewing Results..."
    print results


    print "Exiting Perturbation"
    return results

# sigmaj produces the deviation for the Kpert function from the various statistical
# properties of the system. It takes past distributions and the index you are currently
# on to get the new sigmaj.
def sigmaj(t, index, N, pThetp, wgtp, pThetilde, wgtilde):
    if t == 0:
        return 2 * numpy.var(pThetp[index])
    else:
        sigmabuild = 0
        for i in range(N):
            byte = 0
            for k in range(N):
                bit = numpy.sqrt(wgtp[i]*wgtilde[k]*(pThetilde[k] - pThetp[i])**2)
                byte = byte + bit
            sigmabuild = sigmabuild + byte

# Takes a prior distribution and samples them based on the cumulative version of a new distribution added by the
# user. This means that in general larger parameters are more likely to be sampled.
def CumulativeDistribution(pThetp):
    pThetn = numpy.sort(pThetp)
    ThetBins = set(pThetn)
    ThetVals = numpy.histogram(pThetn, len(ThetBins))
    ThetBins = ThetVals[1]
    ThetUnits = numpy.array(ThetVals[0]).astype('float')
    print 'bins'
    print ThetBins
    print 'vals'
    print ThetUnits
    ThetProps = (1. * ThetUnits)/sum(ThetUnits)
    CumDist = numpy.cumsum(ThetProps)
    randbit =  round(numpy.random.rand(1,1),3)
    ThetInter = numpy.interp(randbit, numpy.insert(CumDist, 0, 0), ThetBins)
    return [ThetInter]


# Takes a distribution function pThetp that represents the distribution of parameters functionalized, the range
# this function is to act over, the number of samples you would like, the type of evnelope you want to use
# and the variance of the Normal distribution should you want it. It produces N number of samples that
# are accepted into the pThetp function.
def RejectionSampling(type = 'Box', var = 1, center = -666.666):

    range = [-20, 20]

    N = 10000

    #def pThetp(x):
    #    if ((x > -1) & (x < 1)):
    #        return 1
    #    else:
    #        return 0

    def pThetp(x):
        return numpy.exp(-numpy.power(x + 3, 2.) / (2 * numpy.power(3, 2.)))

    if (center == -666.666):
        center = numpy.mean(range)

    index = 1
    attempts = []
    success = []

    if (type == 'Box'):

        while index <= N:
            xval = numpy.random.uniform(numpy.amin(range), numpy.amax(range))
            yval = numpy.random.uniform(0,1)
            attempt = [xval, yval]

            attempts = numpy.append(attempts, attempt)

            if yval <= pThetp(xval):
                success = numpy.append(success, xval)
                index = index + 1

    elif (type == 'Normal'):

        g = lambda x: 1/numpy.sqrt(2*3.14159*var)*numpy.exp(-(x-center)**2/(2*var))
        print("min MAX")
        print(numpy.amin(range))
        print(numpy.amax(range))
        bitrange = numpy.linspace(numpy.amin(range), numpy.amax(range), 100)
        print("bitrange")
        print(bitrange)
        print("Normal Dist")
        print (g(bitrange))
        print(".all()?")
        pThetpall = []
        i = 0
        print (len(bitrange))
        while i < len(bitrange):
            bite = (pThetp(bitrange[i]))
            pThetpall.append(bite)
            i = i + 1
        print(pThetpall)
        print("pThetp Dist")
        print (pThetpall)
        print("Scaling Factors?")
        print (pThetpall/g(bitrange))
        ScalingFactor = numpy.amax(pThetpall/g(bitrange))
        print(ScalingFactor)

        while index <= N:

            xval = numpy.random.uniform(numpy.amin(range), numpy.amax(range))
            yval = numpy.random.uniform(0, ScalingFactor * g(xval))
            attempt = [xval, yval]

            attempts = numpy.append(attempts, attempt)

            if yval <= pThetp(xval):
                success = numpy.append(success, xval)
                index = index + 1

    else:

        print "Unrecognized Type: Options Currently available are \"Box\" and \"Normal\"."
        success = -1

    print "results"
    resultsMoo = numpy.histogram(success, bins = 15)

    print resultsMoo

    numpy.savetxt("var2.csv", success, delimiter=",")

    return success


# From Main for Testing
def input_to_df(input_str, N_Species, N_param):
    input_list = map(float, input_str.split(' '))
    theta_index_name = []
    for x in range(0, N_param):
        theta_index_name.append('theta' + str(x))
    Sinit_index_name = []
    for x in range(0, N_Species):
        Sinit_index_name.append('S' + str(x))

    Sinit = pandas.DataFrame(map(int, input_list[0:N_Species]), index = Sinit_index_name, dtype='object')
    print "Initial Species copy numbers used in the simulation algorithm\n", Sinit

    theta = pandas.DataFrame(input_list[N_Species:N_Species + N_param], index = theta_index_name, dtype='object')
    print "Initial theta values\n", theta

    t_param = pandas.DataFrame(map(float, [input_list[N_Species + N_param]]), index = ['delt'], dtype='object')
    t_param.loc['N_time']=map(int, [input_list[N_Species + N_param+1]])
    print "Simulation time parameters\n", t_param
    dict = {'Sinit': Sinit, 'theta': theta, 't_param': t_param}
    df = pandas.concat(dict)
    return df

