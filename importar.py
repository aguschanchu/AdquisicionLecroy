# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 14:42:14 2017

@author: Agustin
"""
from numpy import *
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
from plotly.graph_objs import *
import plotly.graph_objs as go
import plotly.plotly as py
import csv
res = []
for i in range(1,100):
    res.append([])

with open('CaracterizacionPulsos_Corrida2_anchos.txt','r') as minfiles:
    mincsv = csv.reader(minfiles, delimiter=',')
    for row in mincsv:
        res[int(row[1])].append(float(row[0]))
    print(mean(res[25]))
        

with open('anch25.txt', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for x in res[25]:
        wr.writerow([x])

with open('anch40.txt', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for x in res[40]:
        wr.writerow([x])

with open('anch50.txt', 'w') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for x in res[50]:
        wr.writerow([x])