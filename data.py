import pandas as pd
from constants import *
import matplotlib.pyplot as plt


def read_values(index, list):
	list_values = []
	for l in list:
		list_values.append(df_main.iloc[index][l])
	return list_values


df_main = pd.read_csv('Data/merged.txt', sep=",")

# Plot each parameter from graph
for i, var_plot in enumerate(PLOT):
	var_plot_list = PLOT_LIST[i]
	fig = plt.subplot(111)
	for j in range(20): # Plot first 20 values
		fig.plot(var_plot_list, read_values(j, var_plot_list), label=df_main.iloc[j]['GeneName'])
		if(j%1000==0):
			print(j)

	chartBox = fig.get_position()
	fig.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.6, chartBox.height])
	fig.legend(loc='upper center', bbox_to_anchor=(1.45, 1.1), shadow=True, ncol=1)
	plt.title(var_plot)
	# plt.savefig('Plots/' + var_plot+'all')
	plt.show()