import numpy as np
import pandas as pd
from hmmlearn import hmm

df_main = pd.read_csv('Data/cleaned.txt', sep=",")

# X = df_main[['cdRPKM0','cdRPKM1','cdRPKM2']]
# print("Before dropping",len(X))
# X = X.dropna(axis=0, how='any')
# X.to_csv("Data/cleaned.txt", index= False)
# print("After dropping", len(X))

print(df_main.head(5))
X = df_main.values

remodel = hmm.GaussianHMM(n_components=3, covariance_type="full", n_iter=100)
remodel.fit(X)
Z2 = remodel.predict(X)
print(Z2)