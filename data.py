import pandas as pd
from constants import *
import matplotlib.pyplot as plt

def read_values(index,list):
	list_values = []
	for l in list:
		list_values.append(df_main.iloc[index][l])
	return list_values

df_main = pd.read_csv('Data/merged.txt', sep=",")

fig = plt.subplot(111)
for i in range(20):
	fig.plot(TXREADS_LIST, read_values(i,TXREADS_LIST), label = df_main.iloc[i]['GeneName'])

#Plot with legend
chartBox = fig.get_position()
fig.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
fig.legend(loc='upper center', bbox_to_anchor=(1.45, 1.1), shadow=True, ncol=1)
plt.title(TXREADS)
plt.show()
plt.savefig('Plots/'+TXREADS)