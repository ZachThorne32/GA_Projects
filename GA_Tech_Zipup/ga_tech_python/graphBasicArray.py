#For graphing a 1D array
from matplotlib import pyplot as plt
import numpy as np
import arrayInOut as aIO
def graph(arrayText = "output.txt"):
    #graphing our data
    data = aIO.readArray(fileName = "output.txt")
    graphs = plt.figure()
    plt.plot(data, color = 'k')

    #graphing the analytical distribution
    #gaussianSpace = linSpace
    #B = 1.0/ (310 / 3.157746e5)
    #K = 1


    plt.savefig("visualization.png")

