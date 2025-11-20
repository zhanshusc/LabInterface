import matlab.engine
import numpy as np
import pandas as pd


#Adjust file paths accordingly basedon where TAExperiment is, keeping the r at the start
filepath=r"/Users/huishu/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/nachos/src/TAExperiment.m"
eng = matlab.engine.start_matlab()

eng.addpath(filepath, nargout=0)


dat = eng.TAExperiment(r"/Users/huishu/Downloads/LabInterface/pyrazine_30mM_water_MA_s1_05042023")
eng.workspace['dat'] = dat 
arr = eng.eval("dat.TAMean", nargout=1)
np_arr = np.array(arr)
eng.eval("disp(dat.fileName)", nargout=0)
# When you are done
eng.quit()