# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 15:43:38 2017

@author: Agustin
"""

import serial, time
import ast
import visa
rm = visa.ResourceManager('@py')
inst = rm.open_resource("USB0::1155::30016::SPD00002140064::0::INSTR")
arduino = serial.Serial("/dev/ttyACM0", 9600)

outputdir = "/home/iteda/Dropbox/ITeDA/Sistema de control/autotune/codigo 1.0/"
'''
periodoPerturbacion=10 #Periodo en minutos
dutyCyclePerturbacion=0.4 # 0<x<1
esperaInicial = 0 #La idea es no perturbarlo, mientras hace el autotune, no? En minutos
time.sleep(2)

inicio = time.time()
ultimaModificacionPerturbacion = inicio
dutyModificacion=False
voltaje = 0
'''
#Abre el lector serial. El arudino debería strings en el siguiente forato:
#NSensor, Temp, Tiemo
#(donde temp, está en celcius)
#Ej: [2, 27.83, 12.2]
# Se utiliza "v" para comunicarse con fuente de DC
# Se utiliza "e" para enviar un estado del PID con el formato
# [e,tiempo,kp,ki,kd,setpoint,texto]
inst.write("CH1:VOLTage 0")
while True:
    #Guardamos los datos que entrega el arduino via serial
    linea = arduino.readline().decode('utf-8')
    print(linea)
    with open(outputdir+"output.txt",'a') as file:
        file.write(linea+"\n")
    try:
    	linea = ast.literal_eval(linea)
    except:
        linea = ""
    if type(linea) == list:
        #Veamos por casos, que tipo de informacion recibimos
        if type(linea[0]) == int:
        with open(outputdir+"sensor"+str(linea[0])+".txt",'a') as file:
            file.write(str(linea[1])+','+str(linea[2])+"\n")
        elif linea[0] == "v":
            #Cambiamos valor de voltaje, y almacenamos únicamente cuando es alterado
            inst.write("CH1:VOLTage "+str(linea[1]))
            with open(outputdir+"voltajeDC"+".txt",'a') as file:
                file.write(str(linea[1])+','+str(linea[2])+"\n")
        elif linea[0] == "e":
            with open(outputdir+"estado"+".txt",'a') as file:
                lineaw = str(linea[1])
                for i in range(2,7):
                    lineaw=lineaw+','+str(linea[i])
                file.write(lineaw+"\n")
else:
    print(linea)
arduino.close()
