import pandas as pd
from constants import *
import matplotlib.pyplot as plt

max_val = 0 #global variable to keep track of max

def read_values(index, list):
	list_values = []
	for l in list:
		list_values.append(df_main.iloc[index][l])
	local_max = max(list_values)
	global max_val
	if(local_max>max_val):
		max_val = local_max
	return list_values


df_actual = pd.read_csv('Data/merged.txt', sep=",")
df_main = pd.read_csv('Data/cleaned.txt', sep=",")

print(len(df_actual), len(df_main))

# Plot each parameter from graph
for i, var_plot in enumerate(PLOT):
	var_plot_list = PLOT_LIST[i]
	fig = plt.subplot(111)
	for j in range(len(df_main)): # Plot first 20 values
		fig.plot(var_plot_list, read_values(j, var_plot_list))
		if(j%1000==0):
			print(j)
	# chartBox = fig.get_position()
	# fig.set_position([chartBox.x0, chartBox.y0, chartBox.width * 0.6, chartBox.height])
	# fig.legend(loc='upper center', bbox_to_anchor=(1.45, 1.1), shadow=True, ncol=1)
	plt.title(var_plot)
	plt.savefig('Plots/' + var_plot+'all')
	plt.show()