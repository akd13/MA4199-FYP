import numpy as np
import pandas as pd
from hmmlearn import hmm
import warnings
warnings.filterwarnings("ignore")

df_main = pd.read_csv('Data/cleaned.txt', sep=",")

# X = df_main[['cdRPKM0','cdRPKM1','cdRPKM2']]
# print("Before dropping",len(X))
# X = X.dropna(axis=0, how='any')
# X.to_csv("Data/cleaned.txt", index= False)
# print("After dropping", len(X))

print(df_main.head(5))
X = df_main.values
print(X)

df_new = df_main[['cdRPKM0','cdRPKM0']]
print(df_new.head(5))
X1 = df_new.values

remodel = hmm.GaussianHMM(n_components=3,n_iter=100) #consider Gaussian Emissions
remodel.fit(X) #get sequence of states using Viterbi
Z = remodel.predict(X)

# remodel1 = hmm.GaussianHMM(n_components=2, covariance_type="full", n_iter=100)
# remodel1.fit(X1)
# Z1 = remodel.predict(X1)

print("Hidden states are")
for idx,i in enumerate(Z):
	print("State ", idx, " is ",i)
# print(Z1[1:15])

