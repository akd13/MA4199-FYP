# # Import Libraries

# In[40]:


import numpy as np
import pandas as pd
from hmmlearn import hmm
import warnings
from constants import *
import math
import random
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")


# # Helper Methods

def convert_values_to_list(list_val):
	'''
	Given a list X = [1 2 3] ,
	return X = [[1],[2],[3]]
	'''
	X = []
	for i in list_val:
		X.append([i])
	return X


def conversion_list_of_list(X, DIMENSION):
	'''
	Given a list X with values in lists,
	X = [[ 3  2  2], [4 8 10]]
	Convert each value to a list
	Return list of lists, array of lengths of each sequence
	X = [[3] [2] [2] [4] [8] [10]]
	lengths = [3,3]
	'''
	X_new = []
	length = len(X)
	for idx, val_list in enumerate(X):
		Y = []
		for val in val_list:
			Y.append([val])
		X_new.append(Y)

	# Concatenation
	Z = X_new[0]
	for val_list in X_new[1:]:
		Z = np.concatenate([Z, val_list])

	# assign array of lengths for HMM
	lengths = [DIMENSION] * length

	return Z, lengths


# Calculate likelihood for given sequence according to given HMMs and return HMM
def likelihood_sequence(sequence, HMM_array):
	'''
	Given list of K HMMs and sequence = [2 4 6]
	determines likelihood of sequence under all HMM models
	Returns index of HMM which has max likelihood
	'''
	scores = []
	length = [len(sequence)]
	for i, HMM in enumerate(HMM_array):
		calculated_score = HMM.score(sequence, length)
		scores.append(calculated_score)
	idx = scores.index(max(scores))
	return idx


def HMM_model_stats(model):
	'''
	Details of HMM model
	'''
	print("*************************************")
	print("Transition matrix")
	print(model.transmat_)
	print("*************************************")
	print("Means and stds of each hidden state")
	for i in range(model.n_components):
		print("Hidden state {0}".format(i))
		print("mean = ", model.means_[i])
		print("std = ", [math.sqrt(model.covars_[i])])
		print()


def BIC(HMM, X):
	LogLikelihood = model.score(X)
	num_hidden_states = model.n_components
	# D counts transition matrix, emission matrix, sequences estimated (Z), covariance matrix
	D = (num_hidden_states) + 2 * (num_hidden_states ** 2) + len(X) * DIMENSION
	BIC = LogLikelihood - (D / 2) * np.log(len(X))
	return BIC


def BIC_array(HMM_array, X_i):
	BIC_total = 0
	for i in range(len(X_i)):
		model = HMM_array[i]
		X, lengths = conversion_list_of_list(X_i[i])
		LogLikelihood = model.score(X, lengths)
		num_hidden_states = model.n_components
		# D counts transition matrix, emission matrix, sequences estimated (Z), covariance matrix
		D = num_hidden_states + 2 * (num_hidden_states ** 2) + len(X) * DIMENSION
		BIC = LogLikelihood - (D / 2) * np.log(len(X))
		BIC_total += BIC
	return BIC


def plot_BIC(list_k, BIC_score):
	fig = plt.subplot(111)
	plt.plot(list_k, BIC_score, marker='o')
	plt.xlabel('Value of K')
	plt.ylabel('Objective')
	plt.title('BIC')
	plt.show()


# Check increasing sequence
def monotonic_increase(x):
	dx = np.diff(x)
	return np.all(dx >= 0)


# Check decreasing sequence
def monotonic_decrease(x):
	dx = np.diff(x)
	return np.all(dx <= 0)


# Check increasing then decreasing sequence
def increase_decrease(x):
	return x[1] >= x[0] and x[2] <= x[1]


# Check decreasing then increasing sequence
def decrease_increase(x):
	return x[1] <= x[0] and x[2] >= x[1]


# # Load Data and Clean

# In[43]:


# Dataset
df_main = pd.read_csv('Data/cleaned.txt', sep=",")
LENGTH, DIMENSION = df_main.shape
print("Dataset size is", LENGTH)
print("Features are", DIMENSION)
print(df_main.head(5))
X = np.log2(df_main.values)
print("****************************")
print("First 5 log2 values\n", X[:5])

# Generate K HMMs

K = 4
HMM_array = []
X_i = []

for i in range(K):
	X_i.append([])

# Sequences for initial HMM estimation
for i in range(len(X)):
	# Monotone increasing
	if (i % 4 == 0):
		X_i[0].append(list(X[i]))
	elif (i % 4 == 1):
		X_i[1].append(list(X[i]))
	elif (i % 4 == 2):
		X_i[2].append(list(X[i]))
	else:
		X_i[3].append(list(X[i]))


## Train 4 HMMs

for i in range(K):
	model = hmm.GaussianHMM(n_components=3)
	X_temp, lengths = conversion_list_of_list(X_i[i], DIMENSION)
	model.fit(X_temp, lengths)
	HMM_array.append(model)

# # Check likelihood and do assignments

NUM_ITERATIONS = 0
NUM_CLUSTER_PREV = {}
NUM_CLUSTER_NOW = {}

# initialize empty subsets of data
X_i = []

for i in range(K):
	X_i.append([])

while (True):
	# Assign all sequences to HMM models

	print("************ Check likelihood of sequence in HMM  *********")
	NUM_CLUSTER_NOW = {}
	for x in X:
		sequence = convert_values_to_list(x)
		hmm_index = likelihood_sequence(sequence, HMM_array)
		X_i[hmm_index].append(list(x))
		if (hmm_index not in NUM_CLUSTER_NOW):
			NUM_CLUSTER_NOW[hmm_index] = 1
		else:
			NUM_CLUSTER_NOW[hmm_index] += 1
	print("************ Checking likelihood done  *********")


	# Re-estimate parameters for new HMMs
	print("************ Re-estimating HMM *********")
	HMM_array = []
	for i in range(K):
		model = hmm.GaussianHMM(n_components=3)
		X_temp, lengths = conversion_list_of_list(X_i[i], DIMENSION)
		model.fit(X_temp, lengths)
		HMM_array.append(model)
	print("************ Re-estimation done *********")

	# if no reassignments, then break
	if (NUM_CLUSTER_PREV == NUM_CLUSTER_NOW):
		break
	else:
		# initialize empty subsets of data for next iteration
		X_i = []
		for i in range(K):
			X_i.append([])

		NUM_CLUSTER_PREV = NUM_CLUSTER_NOW
		print("Num iterations is:", NUM_ITERATIONS)
		NUM_ITERATIONS += 1
