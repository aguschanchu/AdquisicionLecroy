#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import gc
import psutil
from numpy import *
from plotly.offline import download_plotlyjs, init_notebook_mode,  plot
from plotly.graph_objs import *
import plotly.graph_objs as go
import plotly.plotly as py
import pandas as pd
import scipy
import datetime
import struct
import peakutils
from scipy.integrate import simps
import csv
import datetime
import struct
import sys
import time
import datetime
import math
import argparse
import subprocess
import vxi11
import serial, time
import ast
import visa
import redis

'''
CHANGELOG:
1.0: El codigo para medir, fue definido como una funcion aparte
1.1: El codigo funciona, aunque no entiendo bien por que. Es necesario revisar la documentacion sobre los INR. Que hostia significa INR 0?
1.2: Se agrego un segundo paso al codigo del analisis, teniendo en cuenta los criterios indicados el 19/6 en la reunion con Marie.
Aun persiste el problema de 1.1
1.2.1: Se corrigio un problema en analizador. Ya funciona correctamente.
1.2.2: Integrador ajustado
1.2.3: Se modifico el medidor de anchos, para que sea a mitad de pico
1.3: Se agrego funciones para comunicarse con el sistema de control de Arduino
1.3.1: Nuevo protocolo de comunicacion via Redis
'''
'''

'''

##
#Inicio de readTrc.py
##
def readTrc( fName ):
	"""
		Reads .trc binary files from LeCroy Oscilloscopes.
		Decoding is based on LECROY_2_3 template.
		[More info](http://forums.ni.com/attachments/ni/60/4652/2/LeCroyWaveformTemplate_2_3.pdf)

		Parameters
		-----------
		fName = filename of the .trc file

		Returns
		-----------
		x: array with sample times [s],

		y: array with sample  values [V],

		d: dictionary with metadata


		M. Betz 09/2015
	"""
	with open(fName, "rb") as fid:
		data = fid.read(50).decode()
		wdOffset = data.find('WAVEDESC')

		#------------------------
		# Get binary format / endianess
		#------------------------
		if readX( fid, '?', wdOffset + 32 ):  #16 or 8 bit sample format?
			smplFmt = "int16"
		else:
			smplFmt = "int8"
		if readX( fid, '?', wdOffset + 34 ):  #Big or little endian?
			endi = "<"
		else:
			endi = ">"

		#------------------------
		# Get length of blocks and arrays:
		#------------------------
		lWAVE_DESCRIPTOR = readX( fid, endi+"l", wdOffset + 36 )
		lUSER_TEXT	 = readX( fid, endi+"l", wdOffset + 40 )
		lTRIGTIME_ARRAY  = readX( fid, endi+"l", wdOffset + 48 )
		lRIS_TIME_ARRAY  = readX( fid, endi+"l", wdOffset + 52 )
		lWAVE_ARRAY_1   = readX( fid, endi+"l", wdOffset + 60 )
		lWAVE_ARRAY_2   = readX( fid, endi+"l", wdOffset + 64 )

		d = dict()  #Will store all the extracted Metadata

		#------------------------
		# Get Instrument info
		#------------------------
		d["INSTRUMENT_NAME"]  = readX( fid, "16s",   wdOffset + 76 ).decode().split('\x00')[0]
		d["INSTRUMENT_NUMBER"]= readX( fid, endi+"l", wdOffset + 92 )
		d["TRACE_LABEL"]   = readX( fid, "16s",   wdOffset + 96 ).decode().split('\x00')[0]

		#------------------------
		# Get Waveform info
		#------------------------
		d["WAVE_ARRAY_COUNT"] = readX( fid, endi+"l", wdOffset +116 )
		d["PNTS_PER_SCREEN"]  = readX( fid, endi+"l", wdOffset +120 )
		d["FIRST_VALID_PNT"]  = readX( fid, endi+"l", wdOffset +124 )
		d["LAST_VALID_PNT"]   = readX( fid, endi+"l", wdOffset +128 )
		d["FIRST_POINT"]   = readX( fid, endi+"l", wdOffset +132 )
		d["SPARSING_FACTOR"]  = readX( fid, endi+"l", wdOffset +136 )
		d["SEGMENT_INDEX"]   = readX( fid, endi+"l", wdOffset +140 )
		d["SUBARRAY_COUNT"]   = readX( fid, endi+"l", wdOffset +144 )
		d["SWEEPS_PER_ACQ"]   = readX( fid, endi+"l", wdOffset +148 )
		d["POINTS_PER_PAIR"]  = readX( fid, endi+"h", wdOffset +152 )
		d["PAIR_OFFSET"]   = readX( fid, endi+"h", wdOffset +154 )
		d["VERTICAL_GAIN"]   = readX( fid, endi+"f", wdOffset +156 ) #to get floating values from raw data :
		d["VERTICAL_OFFSET"]  = readX( fid, endi+"f", wdOffset +160 ) #VERTICAL_GAIN * data - VERTICAL_OFFSET
		d["MAX_VALUE"]   = readX( fid, endi+"f", wdOffset +164 )
		d["MIN_VALUE"]   = readX( fid, endi+"f", wdOffset +168 )
		d["NOMINAL_BITS"]   = readX( fid, endi+"h", wdOffset +172 )
		d["NOM_SUBARRAY_COUNT"]= readX( fid, endi+"h",wdOffset +174 )
		d["HORIZ_INTERVAL"]   = readX( fid, endi+"f", wdOffset +176 ) #sampling interval for time domain waveforms
		d["HORIZ_OFFSET"]  = readX( fid, endi+"d", wdOffset +180 ) #trigger offset for the first sweep of the trigger, seconds between the trigger and the first data point
		d["PIXEL_OFFSET"]  = readX( fid, endi+"d", wdOffset +188 )
		d["VERTUNIT"]   = readX( fid, "48s", wdOffset +196 ).decode().split('\x00')[0]
		d["HORUNIT"]   = readX( fid, "48s", wdOffset +244 ).decode().split('\x00')[0]
		d["HORIZ_UNCERTAINTY"]= readX( fid, endi+"f", wdOffset +292 )
		d["TRIGGER_TIME"]  = getTimeStamp( fid, endi, wdOffset +296 )
		d["ACQ_DURATION"]  = readX( fid, endi+"f", wdOffset +312 )
		d["RECORD_TYPE"]   = ["single_sweep","interleaved","histogram","graph","filter_coefficient","complex","extrema","sequence_obsolete","centered_RIS","peak_detect"][ readX( fid, endi+"H", wdOffset +316 ) ]
		d["PROCESSING_DONE"]  = ["no_processing","fir_filter","interpolated","sparsed","autoscaled","no_result","rolling","cumulative"][ readX( fid, endi+"H", wdOffset +318 ) ]
		d["RIS_SWEEPS"]  = readX( fid, endi+"h", wdOffset +322 )
		d["TIMEBASE"]  = ['1_ps/div', '2_ps/div', '5_ps/div', '10_ps/div', '20_ps/rediv', '50_ps/div', '100_ps/div', '200_ps/div', '500_ps/div', '1_ns/div', '2_ns/div', '5_ns/div', '10_ns/div', '20_ns/div', '50_ns/div', '100_ns/div', '200_ns/div', '500_ns/div', '1_us/div', '2_us/div', '5_us/div', '10_us/div', '20_us/div', '50_us/div', '100_us/div', '200_us/div', '500_us/div', '1_ms/div', '2_ms/div', '5_ms/div', '10_ms/div', '20_ms/div', '50_ms/div', '100_ms/div', '200_ms/div', '500_ms/div', '1_s/div', '2_s/div', '5_s/div', '10_s/div', '20_s/div', '50_s/div', '100_s/div', '200_s/div', '500_s/div', '1_ks/div', '2_ks/div', '5_ks/div', 'EXTERNAL'][ readX( fid, endi+"H", wdOffset +324 ) ]
		d["VERT_COUPLING"]  = ['DC_50_Ohms', 'ground', 'DC_1MOhm', 'ground', 'AC,_1MOhm'][ readX( fid, endi+"H", wdOffset +326 ) ]
		d["PROBE_ATT"]  = readX( fid, endi+"f", wdOffset +328 )
		d["FIXED_VERT_GAIN"]  = ['1_uV/div','2_uV/div','5_uV/div','10_uV/div','20_uV/div','50_uV/div','100_uV/div','200_uV/div','500_uV/div','1_mV/div','2_mV/div','5_mV/div','10_mV/div','20_mV/div','50_mV/div','100_mV/div','200_mV/div','500_mV/div','1_V/div','2_V/div','5_V/div','10_V/div','20_V/div','50_V/div','100_V/div','200_V/div','500_V/div','1_kV/div'][ readX( fid, endi+"H", wdOffset +332 ) ]
		d["BANDWIDTH_LIMIT"]  = ['off', 'on'][ readX( fid, endi+"H", wdOffset +334 ) ]
		d["VERTICAL_VERNIER"] = readX( fid, endi+"f", wdOffset +336 )
		d["ACQ_VERT_OFFSET"]  = readX( fid, endi+"f", wdOffset +340 )
		d["WAVE_SOURCE"]   = readX( fid, endi+"H", wdOffset +344 )
		d["USER_TEXT"]  = readX( fid, "{0}s".format(lUSER_TEXT), wdOffset + lWAVE_DESCRIPTOR ).decode().split('\x00')[0]

		#------------------------
		# Get main sample data with the help of numpys .fromfile(
		#------------------------
		fid.seek( wdOffset + lWAVE_DESCRIPTOR + lUSER_TEXT + lTRIGTIME_ARRAY + lRIS_TIME_ARRAY ) #Seek to WAVE_ARRAY_1
		y = fromfile( fid, smplFmt, lWAVE_ARRAY_1 )
		if endi == ">":
			y.byteswap( True )
		y = d["VERTICAL_GAIN"] * y - d["VERTICAL_OFFSET"]
		x = arange(1,len(y)+1)*d["HORIZ_INTERVAL"] + d["HORIZ_OFFSET"]
	return x, y, d

def readX( fid, fmt, adr=None ):
	""" extract a byte / word / float / double from the binary file """
	nBytes = struct.calcsize( fmt )
	if adr is not None:
		fid.seek( adr )
	s = struct.unpack( fmt, fid.read( nBytes ) )
	if(type(s) == tuple):
		return s[0]
	else:
		return s

def getTimeStamp( fid, endi, adr ):
	""" extract a timestamp from the binary file """
	s = readX( fid, endi+"d", adr )
	m = readX( fid, endi+"b" )
	h = readX( fid, endi+"b" )
	D = readX( fid, endi+"b" )
	M = readX( fid, endi+"b" )
	Y = readX( fid, endi+"h" )
	trigTs = datetime.datetime(Y, M, D, h, m, int(s), int((s-int(s))*1e6) )
	return trigTs

##
#Fin de readTrc.py
##

def importar_lecroy(fname):
	'''
	Dado un nombre de archivo .trc, lo importa, y devuelve un array, donde cada una de las entradas, es la traza en formato
	diccionario de Python (claves: tiempo, valores: voltaje)
	'''
	wave = readTrc(fname)
	#nos fijamos cuantas mediciones hay por traza
	mediciones_por_traza = wave[2]["WAVE_ARRAY_COUNT"]//wave[2]["SUBARRAY_COUNT"]
	#llenamos las trazas
	trazas = []
	a = 0
	x=wave[0].tolist()
	y=wave[1].tolist()
	for indice, valor in enumerate(x):
		#si el indice es congruente a 0 modulo numero de trazas, es que empezamos una
		#traza nueva
		if indice % mediciones_por_traza == 0:
			trazas.append(a)
			a = {}
			#paso unidades a mv y ns
			a[float(x[indice])*10**9]=float(y[indice])*1000
		else:
			a[float(x[indice])*10**9]=float(y[indice])*1000
	del x
	del y
	gc.collect()
	return trazas[1:]

def importar_csv_simp(csvf):
	'''
	Dado un archivo csv, de la pinta:
		LECROYWP715Zi-A,57936,Waveform
		Segments,1000,SegmentSize,1002
		Segment,TrigTime,TimeSinceSegment1
		#1,03-May-2017 11:41:40,0
	Lee los parametros (como estan divididas las trazas)
	devuelve un array, donde cada una de las entradas, es la traza en formato
	diccionario de Python (claves: tiempo, valores: voltaje)
	'''

	archivo = csv.reader(csvf, delimiter=',')
	#pasamos el csv a una lista, donde cada objeto es una lista (linea)
	lista = []
	for j in archivo:
		lista.append(j)
	#nos fijamos cuantas mediciones hay por traza
	mediciones_por_traza = int(lista[1][3])
	#quitamos las primeras 3 lineas + #segments +1
	mediciones = lista[3+int(lista[1][1])+1:len(lista)]
	#llenamos las trazas
	trazas = []
	a = 0
	for indice, valor in enumerate(mediciones):
		#si el indice es congruente a 0 modulo numero de trazas, es que empezamos una
		#traza nueva
		if indice % mediciones_por_traza == 0:
			trazas.append(a)
			a = {}
			#paso unidades a mv y ns
			a[float(valor[0])*10**9]=float(valor[1])*1000
		else:
			a[float(valor[0])*10**9]=float(valor[1])*1000
	del lista
	gc.collect()
	return trazas[1:]



def sipm_integrar(traza,graficar=False,anchopico=False,correccion=False):
	'''
	Dada una traza en formato diccionario de Python (claves: tiempo, valores: voltaje)
	toma los primeros 200 puntos para calcular la esperanza y varianza del ruido
	resta el promedio, y busca picos cuyo maximo supera thres
	Luego, busca los limites de integracion para ese pico avanzando en el tiempo
	hasta que V(t) baja por debajo de thres (idem limite inferior)
	Luego, integra numericamente usando simpson
	Devuleve un array con el valor de la integral en cada pico
	Si ademas anchopico=true, devuelve el ancho del pico en ns
	Esta esta adaptada para los valores de la curva del SiPM
	(nada mas peakutils, el resto es igual)
	'''
	#Ruidolimite determina los puntos que se utiliza para analizar E y sigma de Vbaseline.
	#En este caso, es el primer 0.1 (10 por cierto) de la cantidad de puntos totales
	ruidolimite=int(len(traza.keys())*0.05)
	promedio = mean(list(traza.values())[0:ruidolimite])
	desviacion = std(list(traza.values())[0:ruidolimite])
	#Umbral de deteccion de picos, y integracion
	threspeak = 4*desviacion
	thresint = 2*desviacion
	#con el promedio y desviacion, centramos el grafico
	#y lo multiplicamos por -1
	for key in traza:
		traza[key] = -1*(traza[key] - promedio)
	voltajes = list(traza.values())
	tiempos = list(traza.keys())
	#contamos picos:
	#el min_dist es un poco arbitrario, deberia estar relacionado con el ETT
	indices = peakutils.indexes(array(voltajes),min_dist=len(voltajes)//3,thres=threspeak/max(voltajes))
	#filtramos picos que no superen 2*std
	indicesf = []

	for i in indices:
		if voltajes[i]>threspeak:
			indicesf.append(i)

	#para cada indice, buscamos los limites de integracion
	#desplazandonos hacia los costados, hasta bajar mas de 2*std
	#en limites guardamos los limites de cada indice, y el indice
	#en formato [indicemax,indicecotasup,indicecotainf]
	limites = []
	for i in indicesf:
		cotainf=i
		cotasup=i
		try:
			while voltajes[cotasup]>thresint:
				cotasup+=1
			while voltajes[cotainf]>thresint:
				cotainf-=1
			limites.append([i,cotasup,cotainf])
		except:
			if correccion == True:
				#Encontramos una excepcion, ploteo
				graficar = True
				integrales = False
				limites = False
				print("Error al integrar, ECC activo, procedo a devolver False")
				asd = input("presione cualquier tecla para continuar")
			else:
				integrales = False
				limites = False

	#determinamos anchos a mitad de pico
	limitesanchos = []
	for i in indicesf:
		cotainf=i
		cotasup=i
		#calculamos hasta donde hay que bajar (la mitad del pico)
		thresancho=voltajes[i]/2
		try:
			while voltajes[cotasup]>thresancho:
				cotasup+=1
			while voltajes[cotainf]>thresancho:
				cotainf-=1
			limitesanchos.append([i,cotasup,cotainf])
		except:
			if correccion == True:
				#Encontramos una excepcion, ploteo
				graficar = True
				integrales = False
				limites = False
				print("Error al integrar, ECC activo, procedo a devolver False")
				asd = input("presione cualquier tecla para continuar")
			else:
				integrales = False
				limites = False

	#integramos usando el metodo de simpson
	integrales=[]
	if limites != False:
		for i in limites:
			try:
				integrales.append(simps(voltajes[i[2]:i[1]],tiempos[i[2]:i[1]]))
			except:
				if correccion == True:
					#Encontramos una excepcion, ploteo
					graficar = True
					integrales = False
					print("Error al integrar, ECC activo, procedo a devolver False")
					asd = input("presione cualquier tecla para continuar")
				else:
					integrales = False

	#para debuggear, esta bueno graficar:
	if graficar == True:
		trace = go.Scatter(
			x=list(traza.keys()),
			y=voltajes,
			mode='lines',
			name='V(t)',
			)
		trace2 = go.Scatter(
			x=[list(traza.keys())[j] for j in indicesf],
			y=[voltajes[j] for j in indicesf],
			mode='markers',
			marker=dict(
				size=8,
				color='rgb(255,0,0)',
				symbol='cross'
			),
			name='Picos'
			)
		trace3 = go.Scatter(
			x=[list(traza.keys())[ruidolimite]],
			y=[voltajes[ruidolimite]],
			mode='markers',
			marker=dict(
				size=8,
				color='rgb(127,0,255)',
				symbol='cross'
			),
			name='Limite de analisis Vbaseline'
		)

		if limites != False:
			limitest = []
			for i in limites:
				limitest.append(i[1])
				limitest.append(i[2])
			trace4 = go.Scatter(
				x=[list(traza.keys())[j] for j in limitest],
				y=[voltajes[j] for j in limitest],
				mode='markers',
				marker=dict(
					size=8,
					color='rgb(0,255,0)',
					symbol='cross'
				),
				name='Limites de integracion'
			)
			limitesta = []
			for i in limitesanchos:
				limitesta.append(i[1])
				limitesta.append(i[2])
			trace5 = go.Scatter(
				x=[list(traza.keys())[j] for j in limitesta],
				y=[voltajes[j] for j in limitesta],
				mode='markers',
				marker=dict(
					size=8,
					color='rgb(255,128,0)',
					symbol='cross'
				),
				name='Linea de medicion de ancho'
			)

		layout = go.Layout(
		xaxis=dict(
			title='Tiempo (ns)',
				  titlefont=dict(
				family='Courier New, monospace',
				size=18,
				color='#7f7f7f'
			)
		),
		yaxis=dict(
			title='Voltaje (mV)',
				  titlefont=dict(
				family='Courier New, monospace',
				size=18,
				color='#7f7f7f'
			)
		)
		#si queres guardar png, descomenta lo siguiente
		#,width=800*2, height=640*2
		)
		if limites != False:
			data = [trace, trace2, trace3, trace4, trace5]
		else:
			data = [trace, trace2]
		fig = go.Figure(data=data, layout=layout)
		plot(fig)


	if anchopico == True:
		anchos = []
		risetime = []
		if limites != False:
			for i in limitesanchos:
				anchos.append(tiempos[i[1]]-tiempos[i[2]])
			for i in limites:
				risetime.append(tiempos[voltajes.index(max(voltajes))]-tiempos[i[2]])
		else:
			return False, False, False, False
		return integrales, anchos, desviacion, risetime


	return integrales, desviacion

def medir_osc(CH_SCALE,CH_OFFSET,TRIG_LVL,NUM_SEQ,nombre_a_guardar=False,analizar=False):
	'''
	Se conecta al osciloscopio, y mide en el canal 2
	Los argumentos son
	1) Escala vertical [voltaje/div, en string y con unidades]. Ej: "1.0MV"
	2) Offset vertical (misma notacion que 1)
	3) Nivel de trigger (misma notacion que 1)
	4) Numero de secuencias (o trazas) a medir [INT] Ej. 1000
	5) nombre_a_guardar es el nombre del archivo (sin .TRC)
	6) Si analizar es falso, no devuelve nada (solo guarda en disco, de haberlo solicitado)
	a guardar en el disco. Por default, no guarda [STRING o False] Ej. "1001particula"
	Ojo que las tira en workingdir
	Devuelve las mediciones en formato usual
	A continuacion hay una serie de parametros que pueden ser hardcodeados
	'''
	#IP del osc?
	DSO_IP = "192.168.0.2"
	NEW_TIMEOUT = 600*NUM_SEQ
	# Preparo conexion con el DSO
	print('Direccion IP del instrumento: %s\nConectando con instrumento...' %DSO_IP)
	# Asigno la dirección IP.
	InstDSO = vxi11.Instrument(str(DSO_IP))
	# Consulto por el nombre del instrumento. Esto me sirve para ver si estoy
	# conectado al instrumento o no.
	NombreInst = InstDSO.ask("*IDN?")	# Conexión vía LXI
	print('Conectado mediante VXi11 (LXI) al instrumento: %s\n' %NombreInst)
	InstDSO.write("ACAL OFF")	# Auto-calibración interna ON
	InstDSO.write("TIME_DIV 20NS")
	#Trigger tipo edge, del canal 2, sin hold
	InstDSO.write("TRIG_SELECT EDGE, SR, C2, HT, OFF")
	InstDSO.write("TRIG_DELAY -80NS")
	#Triggerea cuando la senal baja del valor de trigger
	InstDSO.write("C2:TRIG_SLOPE NEG")
	InstDSO.write("C2:TRIG_LEVEL " + TRIG_LVL)
	#El valor medio de la señal es 0, de modo, que usamos acoplamiento DC
	InstDSO.write("C2:TRIG_COUPLING DC")
	# Configuración canal 1
	InstDSO.write("C1:TRACE OFF")
	# Configuración canal 2
	InstDSO.write("C2:TRACE ON")
	InstDSO.write("C2:VDIV " + CH_SCALE)
	InstDSO.write("C2:OFST " + CH_OFFSET)
	InstDSO.write("C2:COUPLING D50")
	InstDSO.write("BWL C2,ON")
	# Configuración canal 3
	InstDSO.write("C3:TRACE OFF")
	# Configuración canal 4
	InstDSO.write("C4:TRACE OFF")
	# Configuro la adquisición con sequence
	InstDSO.write("SEQ ON, " + str(NUM_SEQ))
	# Que no almacene trazas en el disco
	InstDSO.write("STST C2, HDD, AUTO, OFF")
	#Le estamos diciendo que almacene toda la informacion
	InstDSO.write("WFSU SP, 0, NP, 0, FP, 0, SN, 0")
	#Formato de trasmision de informacion
	InstDSO.write("CFMT DEF9,WORD,BIN")
	TimeOutVal = InstDSO.timeout
	print('Valor original de TimeOut: %d' %TimeOutVal)
	InstDSO.timeout = NEW_TIMEOUT
	TimeOutVal = InstDSO.timeout
	print('Valor actual de TimeOut: %d' %TimeOutVal)
	##############
	#FIN DE CONFIGURACION OSCILOSCOPIO
	##############
	# Asigno la dirección IP.
	InstDSO = vxi11.Instrument(str(DSO_IP))
	# Consulto por el nombre del instrumento. Esto me sirve para ver si estoy
	# conectado al instrumento o no.
	NombreInst = InstDSO.ask("*IDN?")	# Conexión vía LXI
	print('Dirección IP del instrumento: %s\nConectando con instrumento...' %DSO_IP)
	print(NombreInst)
	INRStatus = 0
	INRStatus = InstDSO.ask("INR?")
	print('INR: %s' %INRStatus)
	try:
		if(INRStatus == 'INR 8193'):
				print('Señal adquirida y lista para ser guardada.\nDSO Listo para ser disparado.')
		elif(INRStatus == 'INR 8192' or INRStatus == 'INR 0' or INRStatus == 'INR 1'):
				print('DSO Listo para ser disparado')
				InstDSO.write("TRMD SINGLE")
				InstDSO.write("ARM")	# Hago un disparo del trigger en modo "single"
				InstDSO.write("WAIT")
				#print(InstDSO.ask("INR?"))
				time.sleep(2)
	except OSError as e:
		print('DSO no disparado: %s %s' %(e, sys.stderr))
	finally:
		print('DSO disparado...')
	print('Esperando a que finalice la adquisición')
	#STB = InstDSO.ask("*STB?")
	#print('STB: %s' %STB)
	INRStatus = False
	#INRStatus = InstDSO.ask("INR?")	# Limpio el registro
	#----------------> TODO: PREGUNTAR POR CODIGOS INR
	#ACASO 8193 O 1 SIGNIFICA QUE YA MIDIO Y ESTA LISTO PARA BAJAR LA TRAZA?
	#Cambiar el while a lo que sea que signifique que haya medido y este listo para entregar signal
	while (INRStatus != 'INR 1' and INRStatus != 'INR 8193' and INRStatus != 'INR 0'):
		try:
			INRStatus = InstDSO.ask("INR?")
		except:
			print(sys.stderr)
		print(INRStatus)
		time.sleep(2)
	print('Adquisición finalizada correctamente')
	print('Obteniendo trazas del oscilo')
	InstDSO.write("C2:WF? ALL")
	try:
		trazas = InstDSO.read_raw()
	except:
		gc.collect()
		return False
	if nombre_a_guardar!=False:
		with open(script_dir+nombre_a_guardar+'.trc','wb') as C1_Trace:
			C1_Trace.write(trazas)
	with open(script_dir+'temp.trc','wb') as C1_Trace:
		C1_Trace.write(trazas)
	time.sleep(2)
	InstDSO.close()
	print('Fin del programa DSO\n')
	if analizar:
		return importar_lecroy(script_dir+'temp.trc')
	else:
		return True



#####
#Fin de definicion de funciones
#####
#CONFIG
script_dir = '/home/iteda/Dropbox/ITeDA/Scripts/Lecroy/resultados/'
temperaturasAMedir = list(range(20,66)) + list(range(5,20))
ktrazasPorTemperatura = 50
#rm = visa.ResourceManager('@py')
#inst = rm.open_resource("USB0::1155::30016::SPD00002140064::0::INSTR")
r = redis.StrictRedis(host='localhost', port=30000, db=0)
#guarda en formato valor,temp
salida=open('minimos.txt','w')
salida2=open('integrales.txt','w')
salida3=open('anchos.txt','w')
salida4=open('risetime.txt','w')
#Cuantas veces queres correr el programa? Tendras 1000*Corridas_globales trazas
guardar_en_disco=True
analizar=False
medir=True
inicio=time.time()
resultados=[]
descartados=0
total=0
#----------------------------------------------------------------------#
# Programa principal
#----------------------------------------------------------------------#
#Por algun motivo aun desconocido, es necesario tirar la primer traza que hagas
#Ver changelog 1.1
if medir == True:
	calentar = medir_osc('2.0MV','7.0MV','-9MV',100)

for temperatura in temperaturasAMedir:
	r.set("setpoint",temperatura)
	#Antes de medir, esperamos que se estabilice la temperatura
	estabilizado=False
	gc.collect()
	while not estabilizado:
		#Guardamos los datos que entrega el arduino via serial
		control = float(r.get("control").decode('UTF-8'))
		print("Temperatura de control: " + str(control))
		print("Temperatura de medicion: " + str(temperatura))
		if abs(control-temperatura) < 0.2:
			estabilizado=True
			#esperamos 15 de estabilizacion en esa temp
			time.sleep(15*60)
		time.sleep(10)

	for k in range(0,ktrazasPorTemperatura):
		print("Iniciando corrida " +str(k)+" setpoint "+str(temperatura))
		corrida = medir_osc('2.0MV','7.0MV','-9.0MV',1000,str(temperatura).zfill(2)+"_"+str(k).zfill(9),False)
		intentos = 0
		while intentos < 10 and corrida == False:
			corrida = medir_osc('2.0MV','7.0MV','-9.0MV',1000,str(temperatura).zfill(2)+"_"+str(k).zfill(9),False)
			intentos+=1

		if analizar == True:
			#Realizamos una primer pasada donde calculamos todo (minimos, anchos, integrales)
		#Luego, en una segunda pasada, descartamos los que se desvian en ksigma
			#de la esperanza
			#Cada resultado corresponde a una matriz con la siguiente sintaxis
			#[desviacion, minimos, [integrales], [anchos],[risetime]]
			for j in range(0, len(corrida)):
				#integracion y anchos
				integrales, anchos, desviacion, risetime = sipm_integrar(corrida[j],False,True,False)
				#nos quedamos unicamente en los casos donde no hubo problemas
				#para integrar. Llenamos los resultados de esta pasada
				resj=[]
				if integrales != False:
					#analisis de minimos
					resj.append(desviacion)
					resj.append(max(corrida[j].values()))
					resj.append(integrales)
					resj.append(anchos)
					resj.append(risetime)
					resultados.append(resj)
				#print(j/len(corrida)*100)

				corrida[j] = ""
				total+=1
		gc.collect()

	#Analizadas todas las trazas, pasamos a descartar
	if analizar == True:
		anchos_primer=[]
		desviaciones_primer=[]
		for i in resultados:
			for j in i[3]:
				anchos_primer.append(j)
			desviaciones_primer.append(i[0])
		umbral_anchos=4*std(anchos_primer)
		promedio_anchos=mean(anchos_primer)
		umbral_desviacion=6*std(desviaciones_primer)
		promedio_desviacion=mean(desviaciones_primer)
		#print(umbral_anchos, promedio_anchos, umbral_desviacion, promedio_desviacion)
		#Con ello, podemos pasar a escribir los archivos
		for i in resultados:
			#si la medicion es muy ruidosa, la descarto
			guardartr = True
			if not (i[0]<promedio_desviacion+umbral_desviacion and i[0]>promedio_desviacion-umbral_desviacion):
				guardartr = False
				descartados+=1
			for j in i[3]:
				#si el ancho es muy grande, lo descarto
				if not (j<promedio_anchos+umbral_anchos and j>promedio_anchos-umbral_anchos):
					guardartr = False
					descartados+=1
			if guardartr == True:
				salida.write(str(i[1])+","+str(temperatura)+","+str(time.time()-inicio)+'\n')
				for j in i[2]:
					salida2.write(str(j)+","+str(temperatura)+","+str(time.time()-inicio)+'\n')
				for j in i[3]:
					salida3.write(str(j)+","+str(temperatura)+","+str(time.time()-inicio)+'\n')
				for j in i[4]:
					salida4.write(str(j)+","+str(temperatura)+","+str(time.time()-inicio)+'\n')


'''
#Ploteamos los histogramas finales
with open(script_dir + '/integrales.txt') as f:
	content = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
a = [log(float(x.strip())) for x in content]
print(len(a))
plot([go.Histogram(x=a,autobinx=True)], filename='integrales.html')
'''



salida.close()
salida2.close()
salida3.close()
