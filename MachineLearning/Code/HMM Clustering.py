
# coding: utf-8

# # Import Libraries

# In[1]:


import numpy as np
import pandas as pd
from hmmlearn import hmm
import warnings
from constants import *
import math
import random
import matplotlib.pyplot as plt
import dill


# # Helper Methods

# In[2]:


def generate_random_sample(X, size):
    '''
    Given a list X, 
    generate random samples of given size
    '''
    Z_temp = random.sample(list(X), size)
    
    #Concatenation
    Z = [Z_temp[0]]
    for val in Z_temp[1:]:
        Z = np.concatenate([Z,[val]])
    
    return Z

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
    X = [[ 3  2  2],[4 8 10]]
    Convert each value to a list
    Return list of lists, array of lengths of each sequence
    X = [[3] [2] [2] [4] [8] [10]]
    ''' 
    X_new = []
    length = len(X)
    for idx, val_list in enumerate(X):
        Y = []
        for val in val_list:
            Y.append([val])
        X_new.append(Y)

    #Concatenation
    Z = X_new[0]
    for val_list in X_new[1:]:
        Z = np.concatenate([Z,val_list])

    # assign array of lengths for HMM
    lengths = [DIMENSION]*length
    
    return Z,lengths

#Calculate likelihood for given sequence according to given HMMs and return HMM
def likelihood_sequence(sequence, HMM_array):
    '''
    Given list of K HMMs and sequence,
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
        print("std = ", [np.sqrt(model.covars_[i])])
        print()

def BIC(model,X,lengths):
    LogLikelihood = model.score(X,lengths)
    num_hidden_states = model.n_components
    # D counts transition matrix (emission estimated by PDF), means = num_hidden_states  
    # covariance matrix = num_hidden_states
    D = num_hidden_states**2-num_hidden_states + 2*num_hidden_states
    BIC = LogLikelihood - (D/2)*np.log(len(X))
    return BIC
    
def BIC_array(HMM_array,X_i):
    BIC_total = 0
    for i in range(len(X_i)):
        model = HMM_array[i]
        X, lengths = conversion_list_of_list(X_i[i],DIMENSION)
        BIC_total+= BIC(model,X,lengths)
    return BIC_total

def likelihood_array(HMM_array,X_i):
    likelihood_total = 1
    for i in range(len(X_i)):
        if(len(X_i[i])>=HMM_array[i].n_components):
            model = HMM_array[i]
            X, lengths = conversion_list_of_list(X_i[i],DIMENSION)
            LogLikelihood = model.score(X,lengths)
            likelihood = LogLikelihood 
            likelihood_total*= likelihood
    return likelihood_total

def plot_BIC(list_k, BIC_score):
    fig = plt.subplot(111)
    plt.plot(list_k, BIC_score, marker='o')  
    plt.xlabel('Value of K')
    plt.ylabel('Objective')
    plt.title('BIC')
    plt.show() 

def print_stats(assignments,length):
    for i in range(length):
        if(i%100==0):
            print(assignments[i])   


# # Load Data and Clean

# In[3]:


df = pd.read_csv('merged.txt', sep=",", na_values=['-'])
df = df.dropna()
df = df[['AccNum','GeneName','cdReads0','cdReads1','cdReads2','cdRPKM0','cdRPKM1','cdRPKM2']]


# # Filter cdReads

# In[4]:


df = df[(df['cdReads0'] >= 10) & (df['cdReads1'] >= 10) & (df['cdReads2'] >= 10)]


# In[5]:


#Dataset
df_main = df[['cdRPKM0','cdRPKM1','cdRPKM2']]
LENGTH,DIMENSION = df_main.shape
print("Dataset size is",LENGTH)
print("Features are", DIMENSION)
print(df_main.head(5))
X = np.log2(df_main.values)
print("****************************")
print("First 5 log2 values\n",X[:5])


# # Arrays with HMM models for 1<=K<=20

# In[6]:


HMM_K_ARRAYS = []
X_i_K_ARRAYS = []


# # Check likelihood and do assignments

# In[7]:


K_values = range(2,21)


# In[8]:


for K in K_values:
    HMM_array = []
    X_i = []
    print("**************** K =", K ,"************************")
    for i in range(K):
        X_i.append([])
        
    NUM_ITERATIONS = 0
    NUM_CLUSTER_PREV = {}
    NUM_CLUSTER_NOW = {}
    
    # Sequences for initial HMM estimation
    # Make K subsets data of LENGTH
    for i in range(LENGTH):
        for j in range(K):
            if(i%K==j):
                X_i[j].append(list(X[i]))
                NUM_CLUSTER_PREV[i] = j
                
    for i in range(K):
        model = hmm.GaussianHMM(n_components=3,covariance_type='diag')
        X_temp, lengths = conversion_list_of_list(X_i[i],DIMENSION)
        model.fit(X_temp, lengths)
        HMM_array.append(model)
    
    likelihood_prev = likelihood_array(HMM_array,X_i)
    print("Likelihood for iteration",NUM_ITERATIONS,"is",likelihood_prev)
    NUM_ITERATIONS+=1
    while (True):
        # Assign all sequences to HMM models

        print("************ Check likelihood of sequence in HMM  *********")
        NUM_CLUSTER_NOW = {}
        for idx,x in enumerate(X):
            sequence = convert_values_to_list(x)
            hmm_index = likelihood_sequence(sequence, HMM_array)
            X_i[hmm_index].append(list(x))
            NUM_CLUSTER_NOW[idx] = hmm_index
        print("************ Checking likelihood done  *********")

        # Re-estimate parameters for new HMMs
        print("************ Re-estimating HMM *********")
        HMM_array_prev = HMM_array
        HMM_array = []
        for i in range(K):
            model = hmm.GaussianHMM(n_components=3,covariance_type='diag')
            if(len(X_i[i])>=model.n_components):
                X_temp, lengths = conversion_list_of_list(X_i[i], DIMENSION)
                model.fit(X_temp, lengths)
                HMM_array.append(model)
            else:
                HMM_array.append(HMM_array_prev[i])
        print("************ Re-estimation done *********")
        likelihood_curr = likelihood_array(HMM_array,X_i)
        print("Likelihood for iteration",NUM_ITERATIONS,"is",likelihood_curr)
        print("*****************************************")

        # if no reassignments, then break
        if ((NUM_CLUSTER_PREV == NUM_CLUSTER_NOW)):
            HMM_K_ARRAYS.append(HMM_array)
            X_i_K_ARRAYS.append(X_i)
            break
        else:
            # initialize empty subsets of data for next iteration
            X_i = []
            for i in range(K):
                X_i.append([])

            NUM_CLUSTER_PREV = NUM_CLUSTER_NOW
            print("Num iterations is:", NUM_ITERATIONS)
            NUM_ITERATIONS += 1
            likelihood_prev = likelihood_curr
    print("**********************************************************\n\n")


# In[17]:


print(len(HMM_K_ARRAYS))
print(len(X_i_K_ARRAYS))


# In[18]:


dill.dump_session('.HMM_GaussianHMM_3points.db')

