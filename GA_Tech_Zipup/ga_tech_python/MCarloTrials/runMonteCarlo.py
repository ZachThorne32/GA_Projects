import sys
import os
sys.path.append('/home/zach/Desktop/ga_tech_python')
import monteCarlo as MC
import graphPosDist as graphing
#temperature
temp = 310
#scale factor/width modifier
scaler = bbb
#Relative duration of trial
dur = aaa
#graphing the result of our trial
graphing.graph(MC.Main(scaler, dur, temp))
