#quantum mechanical dynamics, time dependent schrodinger equation, mean-field (Hartree-
#-Fock) dynamics. All properties of our system contained in 1 electron reduced density matrix.

import numpy as np
import arrayInOut as aIO

class quantumMechanics():
    def __init__(self, n, tStep, h_site, mf1RDM):
        self.TSTEP = tStep
        self.N = n
        self.H = h_site
        self.gamma = mf1RDM
        self.k1Res = None
        self.k2Res = None
        self.k3Res = None
        self.k4Res = None
        #self.trace = np.trace(self.gamma)
        #self.density = np.count_nonzero(self.gamma)/n**2

    def genHermetian(self):
        n = self.N
        nElem = n**2
        matrix = np.random.rand(n, n)+ 1j*np.random.rand(n,n)
        j = 0
        for i in range(n):
            matrix[j, j] = np.real(matrix[j, j])
            j += 1
        matrix = np.triu(matrix)
        return  matrix + np.conj(matrix.T) - np.diag(np.diag(matrix))

    def step(self):
        self.k1()
        self.k2()
        self.k3()
        self.k4()
        self.gamma = np.add(self.gamma, (self.TSTEP / 6.0 * ( np.add(np.add(self.k1Res, self.k2Res * 2.0), np.add(self.k3Res * 2.0, self.k4Res)))))

    def k1(self):
        self.k1Res = commutator(self.H, self.gamma)

    def k2(self):
        self.k2Res = commutator(self.H, np.add(self.gamma, (self.TSTEP / 2.0 * self.k1Res)))

    def k3(self):
        self.k3Res = commutator(self.H, np.add(self.gamma, (self.TSTEP / 2.0 * self.k2Res)))

    def k4(self):
        self.k4Res = commutator(self.H, np.add(self.gamma, (self.TSTEP * self.k3Res)))

def commutator(A, B):
    return -1j*(np.subtract(np.dot(A, B), np.dot(B, A)))

def Main(n, dur, tStep, nPrint, h_site, mf1RDM):
    myQM = quantumMechanics(n, tStep, h_site, mf1RDM)
    t =0
    curTime = np.zeros(dur)
    consistencyCheck = np.zeros([dur, 3])
    data = np.zeros([dur, n+1])
    for i in range(dur):
    	
    	curTime[i] = t
    	consistencyCheck[i, 0] = t
    	data[i, 0] = t
    	totEn = np.real(np.einsum('ij,ji', myQM.H, myQM.gamma))
    	consistencyCheck[i, 1] = (totEn)
    	consistencyCheck[i, 2] = np.real((np.trace(myQM.gamma)))
    	for j in range(n):
    		data[i, 1 + j] = np.real(np.diagonal(myQM.gamma)[j])
    	quantumMechanics.step(myQM)
    	t += myQM.TSTEP
    
    aIO.writeArray(consistencyCheck, 'consistencyCheck.txt')
    aIO.writeArray(data, 'output.txt')
