#!/Anaconda/python
__author__ = 'Keagan Moo'
import numpy
import scipy
import scipy.integrate
import pandas
import csv
import KimScripts as kim

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

# Prior should take in a distribution and select particles from it according to a predetermined distribution
# in our case uniform.
def prior(param_input):
    print "WARNING THIS FUNCTION IS DEPRECATED AND CURRENTLY REPRESENTED BY Initial_particles FROM KIM SCRIPTS"
    # Replace this with whatever your prior distribution ends up being
    Thet = kim.initial_particles(param_input, 1000)
    return Thet

# This version of prior should return a single value as a result of interpreting the probability of a single
# particle being passed in. The probability produced should simply be the sum of the odds of P1 being equal to
# what it is based on the distribution of p1's in the prior distribution and P2 being what it is and P3
# being what it is and so on and so forth. This will all get normalized eventually.
def priorProbDep(param_input, Thet):
    print "priorProb"
    theMany = prior(param_input)
    print "theMany"
    print theMany
    # This strategy is folly. Try again by building a distribution for each paramater and calculating the probability of
    # each parameter originating from that distribution. Then multiply those odds together.

    # Questions: How should we handle rounding? How do you want to change up the prior issue? I could just replace it
    # With the multivariate normal perturbation kernel.

    # Use len(str(#)) to get its number after an initial rounding to say 2 decimal places then round more based on that?
    # It would make more sense to round based on variability though.
    params = numpy.shape(theMany)[1]
    parts = len(theMany[1].loc['theta'])
    paramaterSet = [[0 for x in range(params)] for y in range(parts)]

    for i in range(numpy.shape(theMany)[1]):
        theSubject = theMany[i]
        for j in range(len(paramaterSet)):
            paramaterSet[j][i] = theSubject.loc['theta'][j]

    probProduct = -1
    for p in range(len(paramaterSet)):
        theCount = 0
        itIs = Thet[0].loc['theta'][j]
        paramater = paramaterSet[p]
        for i in range(len(paramater)):
            toMatch = paramater[i]
            if (round(itIs,-2) == round(toMatch, -2)):
                theCount = theCount + 1
        print "The Count of", p
        print theCount
        theProb = (theCount * 1.0)/len(paramater)
        print "The Prob of", p
        print theProb
        if (probProduct == -1):
            probProduct = theProb
        else:
            probProduct = probProduct * theProb
        #print "A Step"
        # probProduct

    print "The Probability"
    print probProduct
    return probProduct

def priorProb(Thet, means, sigmas):
    #print Thet
    paramaters = Thet.loc['theta']
    probProduct = -1
    #print paramaters
    for i in range(len(paramaters)):
        normal = scipy.stats.norm(means[i], sigmas[i])
        prob = scipy.stats.norm(means[i], sigmas[i]).pdf(paramaters[i])

        if (probProduct == -1):
            probProduct = prob
        else:
            probProduct = probProduct * prob


    #print "probProduct"
    #print probProduct
    return probProduct

def normalizeConditional(x):
    total = x.sum()
    if (total == 0):
        return x
    else:
        return x*len(x)/x.sum()


# weightt (Weight Calculator) takes a tolerance counter value,
# a total number of particles, a previous weight vector,
# a theta distribution, the previous theta distribution,
# a perturbation function, and a proposal distribuiton to
# create a new weighting vector for particles.
def weightt(N, prevWeights, pastWeights, newThet, prevThet, i, t, start, sigmas):
    #print "My sins begin to weigh heavily"
    #newWeights = []
    #for i in range(len(newThet)):
    sum = 0
    #print prevThet
    for j in range(N):
        bit = pastWeights[j] * PerturbationProb(newThet[i], prevThet[j], newThet, prevThet, prevWeights, pastWeights, N, t)
        sum = sum + bit
    #print "sum", i
    #print sum
    weight = priorProb(newThet[i], start, sigmas)/sum
    #numpy.append(newWeights, weight)
    #print "The weight falls aside"
    return weight


# Kpert (K perturbation function) stands in for the current perturbation kernel
# it is defined online, but could be changed if needed. It takes the total number
# of particles, a standard deviation and the current and previous parameter distributions
# (theta).

# Change to multi-variate normal (deprecated)
def Kpert(t, N, sigj, pThet, pThetp):
    product = 1
    for j in range(N):
        bit = numpy.exp(-((pThet[j] - pThetp[j])**2)/((2*sigj)**2))/numpy.sqrt(2*3.14159*sigj)
        product = product * bit
    return product

# Multi variate normal kernel without strict adherence to either PNAS or arxiv papers
# Generates particles, not probabilities
def PerturbationKernel(particles, N):
    #print "in kernel"
    ThetaSet = particles.loc['theta']
    #print "rowCount"
    rowCount = ThetaSet.shape[0]
    print "rowCount"
    print rowCount
    colCount = ThetaSet.shape[1]

    means = []
    # Make means of set shape, change how you add means to it.
    index = 0
    #print "While start"
    while index < rowCount:
        subjectParam = ThetaSet.iloc[[index]]
        currentMean = numpy.mean(subjectParam.iloc[0])
        means = numpy.append(means, currentMean)
        index += 1
        #CovarPrep = pandas.DataFrame([CovarPrep,subjectParam])
        #print "Building Covars..."
        #CovarPrep[label] = subjectParam
    print "ThetaSet"
    print ThetaSet
    Covar = numpy.cov(ThetaSet)


    print "means"
    print means
    print "Covar"
    print Covar
    print N

    windex = 0
    results = pandas.DataFrame()
    #print "windex start"
    while windex < N:

        meansl = means.tolist() #[means[0], means[1]]
        Covarl = Covar.tolist() #[[Covar[0][0],Covar[0][1]],[Covar[1][0],Covar[1][1]]]

        result = numpy.random.multivariate_normal(meansl, Covarl)
        results[windex] = result
        windex += 1
    for i in range(N):
        particle = results[i]
        for j in range(len(results)):
            paramater = particle[j]
            if (paramater < 0):
                results[i][j] = results[i][j] * -1

    print "Exiting Perturbation"
    return results

# This function needs to take in both the new and previous sets of particles
# Then using the posteriorFunk of both this tolerance and the previous posterior
# function to integrate over both sets of thetas.
# p = posterior     d = dimensionality
# p(set of particles | Set of real data) the probability of a set of particles existing given the real data
# f(Set of real Data | Set of particles ) the likely hood of a set of real data resulting from certain parameters in particiles.
def PertSum(newThet, prevThet, newW, prevW, N):

    # Left over functionality from the posterior function based version of the perturbation Kernel
    #primaryFunction = prevPosteriorFunk(prevThet, x)*newPosteriorFunk(newThet,x)*(newThet-prevThet)*numpy.transpose((newThet-prevThet))
    #scipy.integrate.dblquad(primaryFunction, 0, 100000,lambda x: 0, lambda x: 100000)
    firstSum = 0
    for i in range(0, N):
        secondSum = 0
        for k in range(0, N):
            difference = newThet[k].loc['theta'] - prevThet[i].loc['theta']
            subSum = 0
            for subSumDex in range(len(difference)):
                subSum = subSum + difference[subSumDex]*difference[subSumDex]
            byte = prevW[i]*newW[k]*subSum
            secondSum = secondSum + byte
        firstSum = firstSum + secondSum

    return firstSum

# This function should take in a newTheta and a previousTheta and return a probability value.
# The probability that such a new theta would result from the old one.

# ERROR: How do we define the single particle in these calculations that is different? Otherwise the sum will always be the same?
# In both the Pert Sum and the prop calculation below we cycle through newThet but it should just be one value. Does PrevThet
# Also just have one value?
def PerturbationProb(newThet, prevThet, newThetWhole, prevThetWhole, newW, prevW, N, t):
    # newThet should be a particle in question
    # prevThet should be a particle to compare it to from a prior set of particles
    newParams = newThet.loc['theta']
    prevParams = prevThet.loc['theta']
    d = len(newParams)
    calculatedSigmat = PertSum(newThetWhole, prevThetWhole, newW, prevW, N)
    #prop = ((2*3.14159)**(-d/2))*((numpy.linalg.det(calculatedSigmat))**(-0.5))*numpy.exp(-0.5*numpy.transpose((newThet-prevThet))*calculatedSigmat^(-1)*(newThet-prevThet))
    matrixStuff = ((newParams-prevParams))
    arrayStuff = numpy.empty(d)
    for index in range(d):
        numpy.append(arrayStuff, matrixStuff[index]*calculatedSigmat**(-1))
    total = 0
    for jndex in range(d):
        total = total + arrayStuff[jndex] * (newParams-prevThet.loc['theta'])[jndex]
    prop = ((2*3.14159)**(-d/2))*((calculatedSigmat)**(-0.5))*numpy.exp(total)
    #if(t == 1):
        #exit(1)
    return prop

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
    ThetProps = (1. * ThetUnits)/sum(ThetUnits)
    CumDist = numpy.cumsum(ThetProps)
    randbit =  round(numpy.random.rand(1,1),3)
    ThetInter = numpy.interp(randbit, numpy.insert(CumDist, 0, 0), ThetBins)
    return [ThetInter]


# Takes a distribution function pThetp that represents the distribution of parameters functionalized, the range
# this function is to act over, the number of samples you would like, the type of evnelope you want to use
# and the variance of the Normal distribution should you want it. It produces N number of samples that
# are accepted into the pThetp function.
def RejectionSampling(type = 'Box', var = 1, center = -123.456):

    range = [-20, 20]

    N = 10000

    #def pThetp(x):
    #    if ((x > -1) & (x < 1)):
    #        return 1
    #    else:
    #        return 0

    def pThetp(x):
        return numpy.exp(-numpy.power(x + 3, 2.) / (2 * numpy.power(3, 2.)))

    if (center == -123.456):
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
        bitrange = numpy.linspace(numpy.amin(range), numpy.amax(range), 100)
        pThetpall = []
        i = 0
        while i < len(bitrange):
            bite = (pThetp(bitrange[i]))
            pThetpall.append(bite)
            i = i + 1

        ScalingFactor = numpy.amax(pThetpall/g(bitrange))
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

def GenerateCorrelate(alpha = 1, beta = -1, gamma = 2, delta = 3, N = 1000):

    x = numpy.random.normal(alpha,gamma, N)
    y = [None]*N
    for i in range(0,len(x)):
        y[i] = numpy.random.normal(beta*x[i], delta, 1)[0]
    primaryDist = pandas.DataFrame({'x':x, 'y':y})

    ofile = open('TestCell.csv', "wb")
    writer = csv.writer(ofile, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
    for i in range(1000):
        z = [x[i], y[i]]
        writer.writerow(z)

    ofile.close()

    return primaryDist


def LinearTransform(S1, S2, dim, N):
    # dim = the number of parameters each particle possesses
    # N = number of particles being input
    xar = ["P1" "P2" "P3" "P4" "P5" "P6" "P7" "P8" "P9" "P10" "P11" "P12"]
    xar = xar[range(dim)]
    Tstar = [None] * N

    # Do we want to just compare each particle to eachother assuming they will be unrelated?
    # Or do we want to compare each particle to every single other particle?

    for i in range(N):
        Tstack = [None] * dim
        for j in range(dim):
            Tbit = S1[i][j]/S2[i][j]
            Tstack[j] = Tbit
        Tstar[i] = Tstack

    # perhaps we get the initial guess for the multiplier this way and assume B = 0 to start with? Eg no translation
    # Either that or we generate purely random numbers for both. But if we want a prior we have to pick one.
    # Then I guess see which system works the best... Again, against a single particle or the whole distribution.
    # Perhaps limit the number of particles to prevent over fitting.

    # The distribution model is just overly simple here. We would just need to calculate the mean of each distribution
    # then the standard deviation. Then the transform would literally be just adding the difference to every particle,
    # and spreading them out from the center by some amount proportional to the new SD.

    # So that would be 2*n numbers to solve for because each parameter we get has both a mean and a standard deviation
    # The transform would then be to multiply each S1 value by its mean modifier and then somehow stratify them?

    transforms = Tstar

    return transforms

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
    #print "Initial Species copy numbers used in the simulation algorithm\n", Sinit

    theta = pandas.DataFrame(input_list[N_Species:N_Species + N_param], index = theta_index_name, dtype='object')
    #print "Initial theta values\n", theta

    t_param = pandas.DataFrame(map(float, [input_list[N_Species + N_param]]), index = ['delt'], dtype='object')
    t_param.loc['N_time']=map(int, [input_list[N_Species + N_param+1]])
    #print "Simulation time parameters\n", t_param
    dict = {'Sinit': Sinit, 'theta': theta, 't_param': t_param}
    df = pandas.concat(dict)
    return df
