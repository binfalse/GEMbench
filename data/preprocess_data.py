import numpy as np
import pandas as pd
import json

ROUND_DIGS = 6
def rounder (num):
  return round (num, ROUND_DIGS)


table = pd.read_csv ("combined-afr-eor-hallmark.csv")

samples = {}
sid = 0
for sample in table["Sample"].unique():
  samples[sample] = sid
  sid = sid + 1

metric_dict = {
  "AFR": 0,
  "EOR": 1,
  "Hallmark": 2,
  "BlandAltman": 3,
  "Jaccard": 4
}
source_dict = {
  "EMTAB-37": 0,
  "GSE2109": 1,
  "HPA": 2,
  "ProteomeNCI60": 3,
  "ProteomePatients": 4,
  "TCGA": 5
}
imeth_dict = {
  "FASTCORE": 0,
  "GIMME": 1,
  "INIT": 2,
  "iMAT": 3,
}

def metric_shorter (metric):
  return metric_dict[metric]

def source_shorter (source):
  return source_dict[source]

def imeth_shorter (imeth):
  return imeth_dict[imeth]

def list_of_samples (subtable, sample_a, sample_b, imeth):
  l = []
  for index, row in subtable.iterrows ():
    tmp = [
            round (row[imeth], ROUND_DIGS),
            samples[row[sample_a]]
          ]
    
    if sample_b:
      tmp.append (samples[row[sample_b]])
    
    l.append (tmp)
  return l

def do_sumstat (metric, source, imeth):
  
  subtable = table[(table.Score==metric) & (table.Dataset==source) & (abs (table[imeth]) != 1000)]
  # subtable = table[(table.Score==metric) & (table.Dataset==source)]
  q1 = rounder (subtable[imeth].quantile (.25))
  median = rounder (subtable[imeth].quantile (.5))
  q3 = rounder (subtable[imeth].quantile (.75))
  
  interQuantileRange = rounder (q3 - q1)
  min_v = rounder (min (subtable[imeth]))
  max_v = rounder (max (subtable[imeth]))
  whiskersMin = rounder (max(min_v, q1 - interQuantileRange * 1.5))
  whiskersMax = rounder (min(max_v, q3 + interQuantileRange * 1.5))
  
  outliers_min = list_of_samples (subtable[rounder (subtable[imeth]) < whiskersMin][["Sample", imeth]], "Sample", None, imeth)
  outliers_max = list_of_samples (subtable[rounder (subtable[imeth]) > whiskersMax][["Sample", imeth]], "Sample", None, imeth)
  
  inftable = table[(table.Score==metric) & (table.Dataset==source) & (abs (table[imeth]) == 1000)]
  inf_p = list_of_samples (inftable[inftable[imeth] == 1000][["Sample", imeth]], "Sample", None, imeth)
  inf_m = list_of_samples (inftable[inftable[imeth] == -1000][["Sample", imeth]], "Sample", None, imeth)
  
  return [ q1,
           median,
           q3,
           interQuantileRange,
           min_v,
           max_v,
           whiskersMin,
           whiskersMax,
           outliers_min,
           outliers_max,
           metric_shorter (metric),
           source_shorter (source),
           imeth_shorter (imeth),
           inf_p,
           inf_m ]




metric = "AFR"
source = "GSE2109"
imeth = "INIT"



sumstat = {}

for metric in ["AFR", "EOR", "Hallmark"]:
  for source in ["HPA", "EMTAB-37", "ProteomeNCI60", "TCGA", "GSE2109", "ProteomePatients"]:
    for imeth in ["GIMME", "iMAT", "INIT", "FASTCORE"]:
      sumstat[str (metric_shorter (metric)) + "_" + str (source_shorter (source)) + "_" + str (imeth_shorter (imeth))] = do_sumstat (metric, source, imeth)




def do_sumstat2 (metric, source, imeth):
  
  subtable = table[(table.Score==metric) & (table.Dataset==source) & (abs (table[imeth]) != 1000)]
  q1 = rounder (subtable[imeth].quantile (.25))
  median = rounder (subtable[imeth].quantile (.5))
  q3 = rounder (subtable[imeth].quantile (.75))
  
  interQuantileRange = rounder (q3 - q1)
  min_v = rounder (min (subtable[imeth]))
  max_v = rounder (max (subtable[imeth]))
  whiskersMin = rounder (max(min_v, q1 - interQuantileRange * 1.5))
  whiskersMax = rounder (min(max_v, q3 + interQuantileRange * 1.5))
  
  outliers_min = list_of_samples (subtable[rounder (subtable[imeth]) < whiskersMin][["Dx1", "Dx2", imeth]], "Dx1", "Dx2", imeth)
  outliers_max = list_of_samples (subtable[rounder (subtable[imeth]) > whiskersMax][["Dx1", "Dx2", imeth]], "Dx1", "Dx2", imeth)
  
  inftable = table[(table.Score==metric) & (table.Dataset==source) & (abs (table[imeth]) == 1000)]
  inf_p = list_of_samples (inftable[inftable[imeth] == 1000][["Dx1", "Dx2", imeth]], "Dx1", "Dx2", imeth)
  inf_m = list_of_samples (inftable[inftable[imeth] == -1000][["Dx1", "Dx2", imeth]], "Dx1", "Dx2", imeth)
  
  return [ q1,
           median,
           q3,
           interQuantileRange,
           min_v,
           max_v,
           whiskersMin,
           whiskersMax,
           outliers_min,
           outliers_max,
           metric_shorter (metric),
           source_shorter (source),
           imeth_shorter (imeth),
           inf_p,
           inf_m ]


table = pd.read_csv ("combined-jaccard-ba.csv")

for metric in ["Jaccard","BlandAltman"]:
  for source in ["HPA", "EMTAB-37", "ProteomeNCI60", "TCGA", "GSE2109", "ProteomePatients"]:
    for imeth in ["GIMME", "iMAT", "INIT", "FASTCORE"]:
      sumstat[str (metric_shorter (metric)) + "_" + str (source_shorter (source)) + "_" + str (imeth_shorter (imeth))] = do_sumstat2 (metric, source, imeth)


















sample_arr = [None] * len (samples)
for i in samples:
  sample_arr[samples[i]] = i
  #print (i)

with open('data.json', 'w') as outfile:
  json.dump({
    "samples": sample_arr,
    "sumstat": sumstat
  }, outfile)
