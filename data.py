import pandas as pd
from constants import *
import matplotlib.pyplot as plt

df_main = pd.read_csv('Data/T0.txt', sep=",")
df_2 = pd.read_csv('Data/T1.txt', sep=",")
df_3 = pd.read_csv('Data/T2.txt', sep=",")

# Merge dataframes

#Merge T1
df_main[TXREADS + '1'] = df_2[TXREADS + '1']
df_main[TXRPKM + '1'] = df_2[TXRPKM + '1']
df_main[CDREADS + '1'] = df_2[CDREADS + '1']
df_main[CDRPKM + '1'] = df_2[CDRPKM + '1']

#Merge T2
df_main[TXREADS + '2'] = df_3[TXREADS + '2']
df_main[TXRPKM + '2'] = df_3[TXRPKM + '2']
df_main[CDREADS + '2'] = df_3[CDREADS + '2']
df_main[CDRPKM + '2'] = df_3[CDRPKM + '2']

df_main.to_csv(path_or_buf='merged.txt', index=False)
print(list(df_main))
print(df_main.iloc[1]['txReads0'])

print(len(df_main))
# plt.plot(['txReads0', 'txReads1', 'txReads2'],[1,2,3])
# plt.show()
