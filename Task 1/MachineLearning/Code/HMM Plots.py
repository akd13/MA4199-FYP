
# coding: utf-8

# In[21]:


from __future__ import with_statement
from __future__ import absolute_import
from io import open
import numpy as np
import pandas as pd
from hmmlearn import hmm
import warnings
from constants import *
import math
import seaborn as sns
import random
import matplotlib.pyplot as plt
warnings.filterwarnings(u"ignore")
import dill
dill.load_session(u'../Weights/HMM_GaussianHMM_3points_RNASeq.db')

list_k = xrange(2,len(X_i_K_ARRAYS)+2)
BIC_scores = []
for idx in xrange(len(X_i_K_ARRAYS)):
    model = HMM_K_ARRAYS[idx]
    X = X_i_K_ARRAYS[idx]
    BIC_scores.append(BIC_array(model,X))


# # In[ ]:
#
#
# for idx,i in enumerate(BIC_scores):
#     print idx+2,i
#
#
# # In[9]:


plot_BIC(list_k,BIC_scores)


# In[ ]:

#
# def plot_cluster(X,count):
#     fig = plt.subplot(111)
#     axes = plt.gca()
#     axes.set_ylim([0,10])
#     var_plot_list = [u'TE0',u'TE1',u'TE2',u'TE3',u'TE4']#['cdRPKM0','cdRPKM1','cdRPKM2','cdRPKM3','cdRPKM4']
#     total=0
#     for i in xrange(len(X)):
#         fig.plot(var_plot_list, X[i])
#         total+=1
#     title = u"HMM "+ unicode(count)+u" : " + unicode(total) + u" points "
#     plt.title(title)
#     plt.savefig('RPKMoutput/Clusters/5points/TE/HMM'+str(count)+'.png')
#     plt.show()
#
#
# # In[ ]:
#
#
# for idx,X in enumerate(X_i_K_ARRAYS[6]):
#     plot_cluster(X,idx+1)
#
#
# # In[ ]:
#
#
# def plot_heatmap(X,idx):
#     plt.figure()
#     sns.heatmap(X,vmin=0, vmax=10)
#     plt.title(u'Heatmap'+unicode(idx))
# #     plt.savefig('RPKMoutput/Clusters/5points/TE/Heatmap'+str(idx+1)+'.png')
#     plt.show()
#
#
# # In[ ]:
#
#
# for idx,X in enumerate(X_i_K_ARRAYS[9]):
#     plot_heatmap(X,idx)
#
#
# # In[ ]:
#
#
# X_model=X_i_K_ARRAYS[0]
#
#
# # In[ ]:
#
#
# df = pd.read_csv(u'RPKMOutput/TE.txt', sep=u" ", na_values=[u'-'])
# df = df.dropna()
# df = df[[u'AccNum',u'GeneName',u'TE0',u'TE1',u'TE2',u'TE3',u'TE4']]
#
#
# # In[ ]:
#
#
# df.head()
#
#
# # In[ ]:
#
#
# # df[['TE0','TE1','TE2','TE3','TE4']] = df[['TE0','TE1','TE2','TE3','TE4']].apply(lambda x: np.log2(x))
#
#
# # In[ ]:
#
#
# for idx,x in enumerate(X_model):
#     genes=[]
#     acc_nums=[]
#     for row in x:
#         temp = (df.loc[(df[u'TE0'] == row[0])& (df[u'TE1']== row[1]) & (df[u'TE2']== row[2]) & (df[u'TE3']== row[3]) & (df[u'TE4']== row[4])])
#         if(not temp.empty):
#             genes.append(temp[u'GeneName'].values[0])
#             acc_nums.append(temp[u'AccNum'].values[0])
#     print len(x),len(genes)
#     with open(u'RPKMoutput/Clusters/5points/TE/GO/Gene'+unicode(idx+1)+u'.txt',u'w') as f:
#         for gene in genes:
#             f.write(u"%s\n" % gene)
#     with open(u'RPKMoutput/Clusters/5points/TE/GO/AccNum'+unicode(idx+1)+u'.txt',u'w') as f:
#         for acc_num in acc_nums:
#             f.write(u"%s\n" % acc_num)
#
