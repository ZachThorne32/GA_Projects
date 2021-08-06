#Monte Carlo code to create configurations for a known distribution
import numpy as np
import random as rand
import arrayInOut as aIO

class MC():
    K = 1
    def __init__(self, scaler, dur, temp):
        self.x = 0
        self.dur = dur
        self.B = 1.0 / (temp / 3.157746e5)
        #self.cStep = cStep
        self.scaler = scaler
        self.pot = self.K / 2 * self.x**2
        self.count = 0
        self.data = np.zeros(int(dur))

    def calcPot(self, xVal):
        potential = self.K / 2 * xVal**2
        return potential

    def addX(self, x):
        self.data[self.count] = x
        self.count += 1

    def compareX(self, newX):
        newPot = self.calcPot(newX)
        if (newPot < self.pot):
            self.addX(newX)
            self.x = newX
            self.pot = newPot
            return
        else:
            P = np.exp(-self.B * (newPot - self.pot))
            if (rand.random() <= P):
                self.addX(newX)
                self.x = newX
                self.pot = newPot
                return
            else:
                self.addX(self.x)

    def stepForward(self):
        newX = self.x + ((rand.random() - .5) * self.scaler)
        self.compareX(newX)

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

    return avgHisto

def Main(scaler, dur, temp):
    myMC = MC(scaler, dur, temp)
    #while data is required, step forward.
    for x in range(int(dur)):
        MC.stepForward(myMC)
    
    aIO.writeArray(calcA(myMC.data, 300), "output.txt")
