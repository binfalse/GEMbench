import numpy as np
import pandas as pd

table = pd.read_csv ("combined-afr-eor-hallmark.csv")

sumstat = {}

metric = "Hallmark"
source = "HPA"
imeth = "GIMME"

q1 = table[(table.Score==metric) & (table.Dataset==source)][imeth].quantile (.25)
median = table[(table.Score==metric) & (table.Dataset==source)][imeth].quantile (.5)
q3 = table[(table.Score==metric) & (table.Dataset==source)][imeth].quantile (.75)


interQuantileRange = q3 - q1;
whiskersMin = max(min (table[(table.Score==metric) & (table.Dataset==source)][imeth]), q1 - interQuantileRange * 1.5);
whiskersMax = min(max (table[(table.Score==metric) & (table.Dataset==source)][imeth]), q3 + interQuantileRange * 1.5);

outliers_min = methvalue["samples"].filter (x => (x.scores[scoreId] < whiskersMin) && vals.includes (x.scores[scoreId]));
outliers_max = methvalue["samples"].filter (x => (x.scores[scoreId] > whiskersMax) && vals.includes (x.scores[scoreId]));
