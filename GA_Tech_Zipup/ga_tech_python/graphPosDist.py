from matplotlib import pyplot as plt
import numpy as np
import arrayInOut as aIO
def graph(arrayText = 'output.txt' ):
#todo: input array from text file, declare B and K explicitly, consider gnuplot graph
    histo = aIO.readArray(fileName = 'output.txt')
    graphs = plt.figure()
    plt.xlim(-.15, .15)
    plt.errorbar(histo[:, 0], histo[:, 1], yerr = histo[:, 2], alpha = .9, errorevery = 2)
    gaussianSpace = np.linspace(-.15, .15, 300)

    B = 1.0 / (310 /3.157746e5)
    K = 1

    analyticalPosDist = np.sqrt(B * K / (2 * np.pi)) * np.exp(-B / 2 * K * gaussianSpace ** 2)

    plt.plot(gaussianSpace, analyticalPosDist, color = 'k')

    plt.savefig('visualization.png')
