import pandas as pd
from constants import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


def read_values(index, list):
	list_values = []
	for l in list:
		list_values.append(df_main.iloc[index][l])
	return list_values

def plot_histogram(list, bins, xlabel, ylabel, title):
	plt.hist(list, bins=bins, facecolor='green')
	plt.xlabel(xlabel)
	plt.ylabel(ylabel)
	plt.title(title)
	plt.show()


df_main = pd.read_csv('Data/cleaned.txt', sep=",")

print("Cleaned dataset size is",len(df_main))

# the histogram of the data
bins = range(0,150,10)
xlabel = 'Gene Expression'
ylabel = 'Number of points'
title = 'Timepoint '
for i in TIMEPOINTS:
	plot_histogram(df_main['cdRPKM'+str(i)],bins,xlabel,ylabel,title+str(i))

# Plot each parameter from graph
for i, var_plot in enumerate(PLOT):
	var_plot_list = PLOT_LIST[i]
	fig = plt.subplot(111)
	for j in range(20): # Plot first 20 values
		print(np.log(var_plot_list))
		fig.plot(var_plot_list, read_values(j, var_plot_list))
		if(j%1000==0):
			print(j)
	# chartBox = fig.get_position()
	# fig.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.6, chartBox.height])
	# fig.legend(loc='upper center', bbox_to_anchor=(1.45, 1.1), shadow=True, ncol=1)
	plt.title(var_plot)
	plt.show()