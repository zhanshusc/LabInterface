import sdtfile
import pandas as pd
import matplotlib.pyplot as plt
sdt = sdtfile.SdtFile('tcspc/irf_261ex_261emt261ex_100uw_2143_nodblexc_05202025_correctcuvetteposition.sdt')
sdt1 = sdtfile.SdtFile('tcspc/tryptamine_0p067odat261ex_350em_200uw_50ns_20mins_1938_may20_2025.sdt')

df = pd.DataFrame(sdt.data[0])
df1 = pd.DataFrame(sdt1.data[0])
print(sdt1.times[0])

plt.plot(sdt1.times[0], df1.iloc[0])

plt.yscale('log')
plt.show()
df1.to_csv("data1.csv",index=False )