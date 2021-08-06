import sys
import os
sys.path.append('/home/zach/Desktop/ga_tech_python/quantumMechanics')
import meanFieldDynamics as QM
import make_hams_multi as hams
import hartreefock

#TODO import some visualization

#duration, N * 20
dur = 10000
#size of matrices (number of electron sites) N + 3 (if odd then add one more)
#print frequency
nPrint = 10
#time steps
tStep = .001

#form hamiltonian
Ndots =2
NL = 2
NR = 2
Vg = 0
U = 0
timp = 1
timplead = 1
Vbias = 0
tleads = 1.0
Full = False
n = Ndots + NL + NR
h_site, V_site = hams.make_ham_multi_imp_anderson_realspace(Ndots, NL, NR, Vg, U, timp, timplead, Vbias, tleads, Full)
#create initial gamma. At n = total electrons is half full
mf1RDM = hartreefock.rhf_calc_hubbard(n, h_site)

Vbias = -.001
h_site, V_site = hams.make_ham_multi_imp_anderson_realspace(Ndots, NL, NR, Vg, U, timp, timplead, Vbias, tleads, Full)

#TODO add visualization step
QM.Main(n, dur, tStep, nPrint, h_site, mf1RDM)

