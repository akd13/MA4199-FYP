from __future__ import with_statement
from __future__ import absolute_import
from google.colab import drive
from io import open
import pandas as pd
import numpy as np
import json
import os
import matplotlib.pyplot as plt

path = sys.argv[0] #/home/akankshitad/Processing/Preprocessing/Data/ATCACG-s_7_1_sequence/
gene = sys.argv[1] #ATCACG-s_7_1_genome
path_to_save = sys.argv[2] #/home/akankshitad/Processing/Preprocessing/Data/Task2/
df = pd.read_csv(path+gene+'.36.bl',sep=u'	',header=None) #'gdrive/My Drive/Task 2/BLFiles/ATCACG-s_7_1_genome.36.bl'

for i in xrange(18,35):
    df_curr=pd.read_csv(path+gene+"."+str(i)+'.bl',sep=u'	',header=None)
    df = df.append(df_curr)

df = df[ (df[10]!= u'rRNA') & (df[11]!=u'tRNA')]

df_known_genes = pd.read_csv('/home/hguo/Documents/annotations/hg19/refFlat/refFlat_240118.txt',sep=u'	',header=None)

for unichr in xrange(22,23):
  position_read = {}
  df_new = df[ df[1]== u'chr'+str(unichr)]
  X = df_new[[8,9]].values
  for x in X:
      if((not np.isnan(x[0])) and (not np.isnan(x[1]))):
          for i in xrange(int(x[0]),int(x[1])+1):
              if(i not in position_read):
                position_read[i]=1
              else:
                position_read[i]+=1
  df_known_chr = df_known_genes[ df_known_genes[2]== u'chr'+str(unichr)]
  X = df_known_chr[[9,10]].values
  known_positions = set()
  for x in X:
    start = x[0].split(u',')[:-1]
    end = x[1].split(u',')[:-1]
    num_pos = len(start)
    for i in xrange(num_pos):
      known_start = int(start[i])
      known_end = int(end[i])
      for j in xrange(known_start, known_end+1):
        if j not in known_positions:
          known_positions.add(j)
  range_chromosome_positions = max(position_read)-min(position_read)
  wiggle_known = [0] * range_chromosome_positions
  wiggle_unknown = [0] * range_chromosome_positions
  for i in xrange(1,range_chromosome_positions+1):
    if(i in position_read and i in known_positions):
      wiggle_known[i]=position_read[i]
    if(i in position_read and i not in known_positions):
      wiggle_unknown[i]=position_read[i]
  fig = plt.figure(figsize=(15,15))
  plt.title(u'chr'+str(unichr))
  plt.xlabel(u'Gene Positions')
  plt.ylabel(u'Number of Reads')
  plt.plot(xrange(min(position_read),max(position_read)), wiggle_known, u'r') # plotting t, a separately 
  plt.plot(xrange(min(position_read),max(position_read)), wiggle_unknown, u'b') # plotting t, b separately 
  plt.savefig(path_to_save+gene+'_chr'+str(unichr)+u'.png', dpi=fig.dpi)
  plt.show()
  with open(u'gdrive/My Drive/Task 2/ATCACG-s_7_1/known_ATCACG-s_7_1_chr'+str(unichr)+u'.txt',u'w') as file:
    file.write(json.dumps(wiggle_known))
  with open(u'gdrive/My Drive/Task 2/ATCACG-s_7_1/unknown_ATCACG-s_7_1_chr'+str(unichr)+u'.txt',u'w') as file:
    file.write(json.dumps(wiggle_unknown))