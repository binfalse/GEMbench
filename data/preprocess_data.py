import numpy as np
import pandas as pd
import json

ROUND_DIGS = 6
def rounder_series (num):
	return round (num, ROUND_DIGS)
def rounder (num):
	if np.isinf (num):
		if num > 0:
			return "I"
		return "i"
	return round (num, ROUND_DIGS)

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

def read_csv (filename):
	df = pd.read_csv (filename)						\
			.replace("Inf", np.inf)					\
			.replace(["-Inf", "#NAME?"], -np.inf)	\
			.replace ("NaN", np.nan)
	for i in imeth_dict:
		df[i] = df[i].astype (float)
	return df

table = read_csv ("combined-afr-eor-hallmark.csv")


samples = {}
sid = 0
for sample in table["Sample"].unique():
  samples[sample] = sid
  sid = sid + 1


def metric_shorter (metric):
  return metric_dict[metric]

def source_shorter (source):
  return source_dict[source]

def imeth_shorter (imeth):
  return imeth_dict[imeth]

def expand_source (source):
  if source == "EMTAB-37":
    return "Cell Line__Microarray"
  if source == "HPA":
    return "Cell Line__RNA Seq"
  if source == "ProteomeNCI60":
    return "Cell Line__MS Proteomics"
  if source == "GSE2109":
    return "Patient Data__Microarray"
  if source == "TCGA":
    return "Patient Data__RNA Seq"
  if source == "ProteomePatients":
    return "Patient Data__MS Proteomics"

def list_of_samples (subtable, sample_a, sample_b, imeth):
  l = []
  for index, row in subtable.iterrows ():
    tmp = [
            rounder (row[imeth]),
            samples[row[sample_a]]
          ]
    
    if sample_b:
      tmp.append (samples[row[sample_b]])
    
    l.append (tmp)
  l.sort(key=lambda x: x[0], reverse=False)
  return l

def do_sumstat (metric, source, imeth):
  
  subtable = table[(table.Score==metric) & (table.Dataset==source) & (abs (table[imeth]) != np.inf)]
  # subtable = table[(table.Score==metric) & (table.Dataset==source)]
  q1 = rounder (subtable[imeth].quantile (.25))
  median = rounder (subtable[imeth].quantile (.5))
  q3 = rounder (subtable[imeth].quantile (.75))
  
  interQuantileRange = rounder (q3 - q1)
  min_v = rounder (min (subtable[imeth]))
  max_v = rounder (max (subtable[imeth]))
  whiskersMin = rounder (max(min_v, q1 - interQuantileRange * 1.5))
  whiskersMax = rounder (min(max_v, q3 + interQuantileRange * 1.5))
  
  outliers_min = list_of_samples (subtable[rounder_series (subtable[imeth]) < whiskersMin][["Sample", imeth]], "Sample", None, imeth)
  outliers_max = list_of_samples (subtable[rounder_series (subtable[imeth]) > whiskersMax][["Sample", imeth]], "Sample", None, imeth)
  
  inftable = table[(table.Score==metric) & (table.Dataset==source) & (abs (table[imeth]) == np.inf)]
  inf_p = list_of_samples (inftable[inftable[imeth] == np.inf][["Sample", imeth]], "Sample", None, imeth)
  inf_m = list_of_samples (inftable[inftable[imeth] == -np.inf][["Sample", imeth]], "Sample", None, imeth)
  
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






def do_sumstat2 (metric, source, imeth):
  
  subtable = table[(table.Score==metric) & (table.Dataset==source) & (abs (table[imeth]) != np.inf)]
  q1 = rounder (subtable[imeth].quantile (.25))
  median = rounder (subtable[imeth].quantile (.5))
  q3 = rounder (subtable[imeth].quantile (.75))
  
  interQuantileRange = rounder (q3 - q1)
  min_v = rounder (min (subtable[imeth]))
  max_v = rounder (max (subtable[imeth]))
  whiskersMin = rounder (max(min_v, q1 - interQuantileRange * 1.5))
  whiskersMax = rounder (min(max_v, q3 + interQuantileRange * 1.5))
  
  outliers_min = list_of_samples (subtable[rounder_series (subtable[imeth]) < whiskersMin][["Dx1", "Dx2", imeth]], "Dx1", "Dx2", imeth)
  outliers_max = list_of_samples (subtable[rounder_series (subtable[imeth]) > whiskersMax][["Dx1", "Dx2", imeth]], "Dx1", "Dx2", imeth)
  
  inftable = table[(table.Score==metric) & (table.Dataset==source) & (abs (table[imeth]) == np.inf)]
  inf_p = list_of_samples (inftable[inftable[imeth] == np.inf][["Dx1", "Dx2", imeth]], "Dx1", "Dx2", imeth)
  inf_m = list_of_samples (inftable[inftable[imeth] == -np.inf][["Dx1", "Dx2", imeth]], "Dx1", "Dx2", imeth)
  
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
  boxplots = {}
  minY = None
  maxY = None
  infinites_p = False
  infinites_m = False
  xdomain = []
  for source in ["HPA", "EMTAB-37", "ProteomeNCI60", "TCGA", "GSE2109", "ProteomePatients"]:
    for imeth in ["GIMME", "iMAT", "INIT", "FASTCORE"]:
      curstat = do_sumstat (metric, source, imeth)
      if minY is None:
        minY = curstat[4]
        maxY = curstat[5]
        infinites_p = len (curstat[13]) > 0
        infinites_m = len (curstat[14]) > 0
      else:
        minY = min (minY, curstat[4])
        maxY = max (maxY, curstat[5])
        infinites_p = infinites_p or len (curstat[13]) > 0
        infinites_m = infinites_m or len (curstat[14]) > 0
      xdomain.append (metric + "_" + expand_source (source) + "_" + imeth)
      boxplots[str (metric_shorter (metric)) + "_" + str (source_shorter (source)) + "_" + str (imeth_shorter (imeth))] = curstat
  sumstat[metric] = {
    "boxplots": boxplots,
    "minY": minY,
    "maxY": maxY,
    "infinites_p": infinites_p,
    "infinites_m": infinites_m,
    "xdomain": xdomain
  }




table = read_csv ("combined-jaccard-ba.csv")

for metric in ["Jaccard","BlandAltman"]:
  boxplots = {}
  minY = None
  maxY = None
  infinites_p = False
  infinites_m = False
  xdomain = []
  for source in ["HPA", "EMTAB-37", "ProteomeNCI60", "TCGA", "GSE2109", "ProteomePatients"]:
    for imeth in ["GIMME", "iMAT", "INIT", "FASTCORE"]:
      curstat = do_sumstat2 (metric, source, imeth)
      if minY is None:
        minY = curstat[4]
        maxY = curstat[5]
        infinites_p = len (curstat[13]) > 0
        infinites_m = len (curstat[14]) > 0
      else:
        minY = min (minY, curstat[4])
        maxY = max (maxY, curstat[5])
        infinites_p = infinites_p or len (curstat[13]) > 0
        infinites_m = infinites_p or len (curstat[14]) > 0
      xdomain.append (metric + "_" + expand_source (source) + "_" + imeth)
      boxplots[str (metric_shorter (metric)) + "_" + str (source_shorter (source)) + "_" + str (imeth_shorter (imeth))] = curstat
  sumstat[metric] = {
    "boxplots": boxplots,
    "minY": minY,
    "maxY": maxY,
    "infinites_p": infinites_p,
    "infinites_m": infinites_m,
    "xdomain": xdomain
  }
      # sumstat[str (metric_shorter (metric)) + "_" + str (source_shorter (source)) + "_" + str (imeth_shorter (imeth))] = do_sumstat2 (metric, source, imeth)


















sample_arr = [None] * len (samples)
for i in samples:
  sample_arr[samples[i]] = i
  #print (i)

with open('data.json', 'w') as outfile:
  json.dump({
    "samples": sample_arr,
    "sumstat": sumstat
  }, outfile)
