import matplotlib.pyplot as plt
import matlab.engine
import numpy as np
import pandas as pd
import seaborn as sns



#Adjust file paths accordingly basedon where TAExperiment is, keeping the r at the start
filepath=r"./TAExperiment.m"
eng = matlab.engine.start_matlab()

eng.addpath(filepath, nargout=0)


dat = eng.TAExperiment(r"./pyrazine_30mM_water_MA_s1_05042023")
eng.workspace['dat'] = dat 
arr = eng.eval("dat.TAMean", nargout=1)
times = eng.eval("dat.times", nargout=1)
wavelengths = eng.eval("dat.wavelengths",nargout=1)
np_arr = np.array(arr)
np_times =np.array(times).flatten()
wavelengths = np.array(wavelengths).flatten()
df = pd.DataFrame(np_arr, index=times, columns=wavelengths)
plt.figure(figsize=(12, 8))
sns.heatmap(df, xticklabels=50, yticklabels=40, cmap='viridis')
plt.show()
eng.eval("disp(dat.fileName)", nargout=0)
# When you are done
eng.quit()