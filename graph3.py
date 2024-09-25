import time
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
import os

#Uso de Latex en los ejes
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams['font.size'] = 14
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "serif"

#Leer Archivo de N = 4
with open("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p3/N4_magnetization_vs_temperature.dat", "r") as f1:
    lines1 = f1.readlines()
    T = [float(line.split()[0]) for line in lines1]
    M_N4 = [float(line.split()[1]) for line in lines1]

#Leer Archivo de N = 8
with open("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p3/N8_magnetization_vs_temperature.dat", "r") as f2:
    lines2 = f2.readlines()
    M_N8 = [float(line.split()[1]) for line in lines2]
    
#Leer Archivo de N = 16
with open("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p3/N16_magnetization_vs_temperature.dat", "r") as f3:
    lines3 = f3.readlines()
    M_N16 = [float(line.split()[1]) for line in lines3]
    
#Leer Archivo de N = 32
with open("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p3/N32_magnetization_vs_temperature.dat", "r") as f4:
    lines4 = f4.readlines()
    M_N32 = [float(line.split()[1]) for line in lines4]

#Grafico de Magnetizacion vs Temperatura
plt.xlabel("$T$")
plt.ylabel(r"$\langle M \rangle$")
plt.plot(T, M_N4, label="N=4")
plt.plot(T, M_N8, label="N=8")
plt.plot(T, M_N16, label="N=16")
plt.plot(T, M_N32, label="N=32")
plt.grid(alpha =.6, linestyle ='--')
plt.legend()
plt.show()
plt.savefig("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/latex/p3/TM.pdf")