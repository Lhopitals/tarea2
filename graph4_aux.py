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

#Leer Archivo
with open("p4/N6_binder_vs_temperature.dat", "r") as f1:
    lines1 = f1.readlines()
    T = [float(line.split()[0]) for line in lines1]
    binder = [float(line.split()[1]) for line in lines1]
    
with open("p4/N16_binder_vs_temperature.dat", "r") as f2:
    lines2 = f2.readlines()
    T = [float(line.split()[0]) for line in lines2]
    binder2 = [float(line.split()[1]) for line in lines2]
    
with open("p4/N32_binder_vs_temperature.dat", "r") as f3:
    lines3 = f3.readlines()
    T = [float(line.split()[0]) for line in lines3]
    binder3 = [float(line.split()[1]) for line in lines3]

with open("p4/N64_binder_vs_temperature.dat", "r") as f4:
    lines4 = f4.readlines()
    T = [float(line.split()[0]) for line in lines4]
    binder4 = [float(line.split()[1]) for line in lines4]





#Grafico de Magnetizacion vs Temperatura
plt.xlabel("$T$")
plt.ylabel(r"$B(T)$")
# plt.plot(T, binder, label="N=6")
# plt.plot(T, binder2, label="N=16")
# plt.plot(T, binder3, label="N=32")
plt.plot(T, binder4, label="N=64")
# plt.axvline(x=2/(np.log(1+np.sqrt(2))), color='black', linestyle='--', label="$T_c$")
# plt.plot(T, binder2, label="N=4")
# plt.plot(T, M_N16, label="N=16")
# plt.plot(T, M_N32, label="N=32")
plt.grid(alpha =.6, linestyle ='--')
plt.legend()
plt.show()
# plt.savefig("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/latex/p2/N2_binder.pdf")