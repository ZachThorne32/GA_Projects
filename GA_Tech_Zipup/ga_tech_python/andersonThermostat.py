#!/usr/bin/env python
# coding: utf-8

# In[1]:


from matplotlib import pyplot as plt
import random as rand
import numpy as np
import arrayInOut as aIO

#Importing matplotlib.pyplot for graphing purposes, numpy for math functions and array manipulation
#and importing random for gaussian distribution.

#Setting e and pi to themselves.
e = np.e
pi = np.pi
#The goal of V2 is to change my data collection structure to numpy arrays (Done), separate the tempSample call from
#collection of important data as to allow for more logical data collection (Done), more explicitly define min and max values
#for the histogram bounds in preparation for error bars (Done), and utilize more numpy packages if possible.
#Next step is calculating the standard of the mean in order to make error bars.


# In[2]:


def forces(x, K):
    #Calculates and returns the force at position x.
    return -K * x


# In[3]:


class MD:
    def __init__(self, k, m, x, v, timeStep, nTemp, T):
        #Initializes properties of the particle and environment parameters. Also calls tempSample in order
        #to start with action
        KB = 1.38e-23
        #B = 1 / (T * KB)
        global B
        global K
        K = k
        global M
        M = m
        B = 1.0 / (T / 3.157746e5) #Changed beta calculations to account for units
        self.f = forces(x, K)
        self.x = x
        self.v = v
        self.mu = 0
        self.sigma = 1 / (np.sqrt(B * M))
        tempSample(self)
        self.timeStep = timeStep
        self.xBackStep = x - v * timeStep
        self.t = 0
        self.nTemp = nTemp


# In[4]:


def positions(myMD):
    #Calculates and returns the new position while updating the backstep (previous) position.
    temp = myMD.x
    newX = myMD.x + myMD.v * myMD.timeStep
    myMD.xBackStep = temp
    myMD.t += myMD.timeStep
    return newX


# In[5]:


def velocities(myMD):
    #Calculates and returns the velocity a half time step ahead of the current position. This will be called twice
    #for each time step. the position call should be using the average velocity for the time step but the velocity
    #still needs to eventually be fully updated, hence the second call.
    newV = myMD.v + 1 / (M * 2) * myMD.f * myMD.timeStep
    return newV


# In[6]:


def tempSample(myMD):
    #Randomly set velocity within a gaussian distribution with mu = 0, standard deviation = 1/sqrt(Beta*Mass)
    #This occurs at variable frequencies set by the nTemp variable
    #The reason we only do this every once in awhile is to have a sample of decorrolated points not corrolated points
    myMD.v = rand.gauss(myMD.mu, myMD.sigma)


# In[7]:


def calcA(array, iBins):
    nChunks = 5
    nBins = iBins
    splitArray = np.split(array, nChunks)
    histo = np.zeros( [nBins, nChunks+1])
    #300 rows of 6 elements, all zero
    for i in range(nChunks): #Runs 0-(bins-1), total of "bins" times.
        if( i == 0 ):
            #First time through create binEdges and calculate midpoints, then fill the histogram with midpoints.
            binEdges = np.histogram( splitArray[i], bins=nBins, range=( array.min(),
            array.max() ), density=True )[1] #assigns the second return ([1]) to binEdges
            for j in range( nBins ): #for 0-(nBins-1)
                histo[j,0] = ( binEdges[j] + binEdges[j+1] ) / 2.0
                #x axis in the histogram gets set to midpoints of the histogram

        #calculate normalized histogram for each chunk, note that all histograms have the same range
        histo[:,i+1]   = np.histogram( splitArray[i], bins=nBins, range=( array.min(), array.max() ), density=True )[0]
    avgHisto = np.zeros( [nBins,3] )
    avgHisto[:,0] = np.copy( histo[:,0] )
    avgHisto[:,1] = np.mean( histo[:,1:], axis=1 )
    avgHisto[:,2] = np.std( histo[:,1:], axis=1 ) / np.sqrt(nChunks-1)

    #print out P(x)
    #return avgHisto[:,0], avgHisto[:, 1], avgHisto[:, 2];
    return avgHisto


# In[8]:


def Main(K, M, xi, vi, Step, nTemp, T, durr):
    #

    myMD = MD(K, M, xi, vi, Step, nTemp, T)

    sampCount = 0

    # on the order of 5,000,000 points for ideal gaussian histogram. Using less to save .
    #Run: .5 mil ~ 6s 1 mil ~ 10s, 5 mil ~ 23s, 50 mil ~ 3 min 17s, 100 mil no line graphs ~ 5:45
    durr
    dataFreq = 10

    Axis = np.linspace(0, 1000, durr)
    gaussianSpace = np.linspace(-.15, .15, 300)
    #posAxis = np.zeros()
    #velAxis = np.zeros()
    distSpace = round(durr/dataFreq)
    posDistribution = np.zeros(distSpace)
    velDistribution = np.zeros(distSpace)
    analyticalPosDist = np.sqrt(B * K / (2 * pi)) * np.exp(-B / 2 * K * gaussianSpace ** 2)
    analyticalVelDist = np.sqrt(B * M / (2 * pi)) * np.exp(-B / 2 * M * gaussianSpace ** 2)

    myMD.f = forces(myMD.x, K)

    for i in range(durr):
        #posAxis[i] = myMD.x
        #velAxis[i] = myMD.v

        myMD.v = velocities(myMD)
        myMD.x = positions(myMD)
        myMD.f = forces(myMD.x, K)
        myMD.v = velocities(myMD)
        if (i % nTemp == 0):
            tempSample(myMD)
        if (i % dataFreq == 0):
            posDistribution[sampCount] = myMD.x
            velDistribution[sampCount] = myMD.v
            sampCount += 1

    #posMin = posAxis.min()
    #posMax = posAxis.max()
    #velMin = velAxis.min()
    #velMax = velAxis.max()

    #posGraph = plt.figure()
    #posGraph.set_figwidth(15)
    #posGraph.set_figheight(4)
    #plt.plot(Axis, posAxis, color = 'g')
    #plt.xlim(0, 400)
    #plt.ylim(-1.1, 1.1)
    #plt.xlabel(' (s)')
    #plt.ylabel('Position')
    #plt.title('Positions Over ')
    #plt.show()

    #velGraph = plt.figure()
    #velGraph.set_figwidth(15)
    #velGraph.set_figheight(4)
    #plt.plot(Axis, velAxis, color = 'c')
    #plt.xlim(0, 400)
    #plt.ylim(-1.1, 1.1)
    #plt.xlabel(' (s)')
    #plt.ylabel('Velocity')
    #plt.title('Velocity Over ')
    #plt.show()
    Ahisto = calcA(posDistribution, len(gaussianSpace))

    #velDist = plt.figure()
    #velDist.set_figwidth(15)
    #velDist.set_figheight(6)
    #plt.xlim(-.15, .15)
    #plt.hist(velAxis, bins = 300, density = True, color = 'm', alpha=0.75)
    #plt.plot(gaussianSpace, analyticalVelDist, color = 'b')
    #plt.xlabel('Velocities')
    #plt.ylabel('Frequency')
    #plt.title('Velocity Distribution')
    #plt.show()


    #calcA(velAxis)
#Models, samples, and graphs a classic moledular dynanamics simulation of harmonic motion
    print('Simulation complete; data available!')
    aIO.writeArray(Ahisto, "output.txt")
#    return Ahisto
#ToDo: write histo to a textfile

#def graph(histo):
#    Ahisto = histo
#    posDist = plt.figure()
#    posDist.set_figwidth(16)
#    posDist.set_figheight(6)
#    plt.xlim(-.15, .15)
    #histoSpace = np.linspace(-.15, .15, len(analyticalPosDist))
    #plt.errorbar(histoSpace, posAxis, yerr = Ahisto, alpha = 0.75, color = 'b', linestyle = '-')
#    plt.errorbar(Ahisto[:, 0], Ahisto[:, 1] , yerr = Ahisto[:, 2], alpha = 0.90, color = 'b', errorevery = 1)
    #plt.plot(time, posAxis, color='c')
#    gaussianSpace = np.linspace(-.15, .15, 300)
#    analyticalPosDist = np.sqrt(B * K / (2 * pi)) * np.exp(-B / 2 * K * gaussianSpace ** 2)
#    plt.plot(gaussianSpace, analyticalPosDist, color = 'k')


