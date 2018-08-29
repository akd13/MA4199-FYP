import pandas as pd
from constants import *
import matplotlib.pyplot as plt

def read_values(index,list):
	list_values = []
	for l in list:
		list_values.append(df_main.iloc[index][l])
	return list_values

df_main = pd.read_csv('Data/merged.txt', sep=",")

for i in range(100):
	plt.plot(TXREADS_LIST, read_values(i,TXREADS_LIST), label = df_main.iloc[i]['GeneName'])

plt.legend(loc = 'upper right')
plt.show()
