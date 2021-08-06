import sys
import os
sys.path.append('/home/zach/Desktop/ga_tech_python')
import andersonThermostat as anderTherm
import graphPosDist as graphing
#Force constant
k = 1
#Mass
m = 1

#Initial position
xi = 0
#Initial velocity
vi = 1
#Time step
tStep = .1
#Frequency of velocity change (interaction with surrounding environment)
nTemp = 250
#Temperature
temp = 310
#Duration of trial
dur = 2000000
#anderTherm.Main(k, m, xi, vi, tStep, nTemp, temp, dur)
graphing.graph(anderTherm.Main(k, m, xi, vi, tStep, nTemp, temp, dur))
