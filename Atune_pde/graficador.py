#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 24 17:21:29 2017

@author: agus
"""

from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
from plotly.graph_objs import *
import plotly.graph_objs as go
import plotly.plotly as py
import time
import gc

while True:
	try:
		dire=''
		res = []
		for i in range(0,9):
			res.append([[],[]])
			with open(dire+"sensor"+str(i)+".txt",'r') as filem:
				for line in filem.readlines():
					line = line.split(',')
					res[i][0].append(float(line[1])/60)
					if i==6:
						res[i][1].append(float(line[0])/255)
					else:
						res[i][1].append(float(line[0]))

		trace1 = go.Scatter(
		    x=res[0][0],
		    y=res[0][1],
		    name="Sensor 0"
		)
		trace2 = go.Scatter(
		    x=res[1][0],
		    y=res[1][1],
		    name="Sensor 1"
		)
		trace3 = go.Scatter(
		    x=res[2][0],
		    y=res[2][1],
		    name="Sensor 2"
		)
		trace4 = go.Scatter(
		    x=res[3][0],
		    y=res[3][1],
		    name="Sensor 3"
		)

		trace5 = go.Scatter(
		    x=res[4][0],
		    y=res[4][1],
		    name="Sensor 4"
		)

		trace6 = go.Scatter(
		    x=res[5][0],
		    y=res[5][1],
		    name="Temperatura ambiente"
		)

		trace7 = go.Scatter(
		    x=res[6][0],
		    y=res[6][1],
		    name="Potencia normalizada"
		)

		trace8 = go.Scatter(
			x=res[7][0],
			y=res[7][1],
			name="Sensor 0 promediado (Control)"
		)

		trace9 = go.Scatter(
			x=res[8][0],
			y=res[8][1],
			name="Potencia resistencia"
		)

		layout = go.Layout(width=800*2, height=640*2,
			xaxis=dict(
					title='Tiempo (minutos)',
						  titlefont=dict(
						family='Courier New, monospace',
						size=18,
						color='#7f7f7f'
					)
				),
				yaxis=dict(
					title='Temperatura (C)',
						  titlefont=dict(
						family='Courier New, monospace',
						size=18,
						color='#7f7f7f'
					))
			)
		data = [trace1,trace2,trace3,trace4,trace5,trace6,trace7,trace8,trace9]
		fig = go.Figure(data=data, layout=layout)
		plot(fig, auto_open=False)
		print("Graficando")
		time.sleep(120)
		gc.collect()
	except:
		pass
