import sdtfile
import pandas as pd
import matplotlib.pyplot as plt
sdt = sdtfile.SdtFile('irf_261ex_261emt261ex_100uw_2143_nodblexc_05202025_correctcuvetteposition (1).sdt')
sdt1 = sdtfile.SdtFile('tryptamine_0p067odat261ex_350em_200uw_50ns_20mins_1938_may20_2025.sdt')

df = pd.DataFrame(sdt.data[0])
df1 = pd.DataFrame 
print(sdt.measure_info[0])

plt.plot(df.iloc[0])


plt.yscale('log')
plt.show()
#df.to_csv("data1.csv",index=False )