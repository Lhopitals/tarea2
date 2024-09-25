import numpy as np
import matplotlib.pyplot as plt

# Uso de Latex en los ejes
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams['font.size'] = 14
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "serif"

# Leer archivo T = 2
with open("p2/b/N4_magnetizacion_distribucion_GrayFlipT2.dat", "r") as f1:
    lines1 = f1.readlines()
    M = [float(line.split()[0]) for line in lines1]  # Magnetización (primera columna)
    Pi_m = [float(line.split()[1]) for line in lines1]  # Probabilidad (segunda columna)

#Leer archivo para T = 5
with open("p2/b/N4_magnetizacion_distribucion_GrayFlipT5.dat", "r") as f2:
    lines2 = f2.readlines()
    M2 = [float(line.split()[0]) for line in lines2]  # Magnetización (primera columna)
    Pi_m2 = [float(line.split()[1]) for line in lines2]  # Probabilidad (segunda columna)

# #Leer archivo T = 2
# with open("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p2/b/N4_E_GrayFlipT2.dat", "r") as f1:
#     lines1 = f1.readlines()
#     E1 = [float(line.split()[0]) for line in lines1]  # Magnetización (primera columna)
#     count1 = [float(line.split()[1]) for line in lines1]  # Probabilidad (segunda columna)
    
# #Leer archivo T = 5
# with open("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p2/b/N4_E_GrayFlipT5.dat", "r") as f2:
#     lines2 = f2.readlines()
#     E2 = [float(line.split()[0]) for line in lines2]  # Magnetización (primera columna)
#     count2 = [float(line.split()[1]) for line in lines2]  # Probabilidad (segunda columna)

# # Leer archivo T = 2
# with open("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p2/b/N4_N(E,M)_GrayFlipT2.dat/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p2/b/N4_energy_GrayFlipT2.dat", "r") as f1:
#     lines1 = f1.readlines()
#     E1 = [float(line.split()[0]) for line in lines1]  # Magnetización (primera columna)
#     count1 = [float(line.split()[1]) for line in lines1]  # Probabilidad (segunda columna)


# Crear histograma ponderado por las probabilidades
plt.xlabel("$M$")
plt.ylabel("$\pi_M$")
plt.hist(M, bins=5, weights=Pi_m, edgecolor='black', alpha=0.7, density=True, histtype="step", label='T = 2', linestyle="--")  # Histograma ponderado por Pi_m
plt.hist(M2, bins=5, weights=Pi_m2, edgecolor='black', alpha=0.7, density=True, histtype="step", label='T = 5')  # Histograma ponderado por Pi_m
# plt.xlabel("$E$")
# plt.xlim(-9, 9)
# plt.xticks([-8, 0 ,8])  # Define los ticks en el eje x
# plt.hist(E1, bins=16, weights=count1, edgecolor='black', alpha=0.7, density=True, histtype="step", label='T = 2', linestyle="--")  # Histograma ponderado por Pi_m
# plt.hist(E2, bins=16, weights=count2, edgecolor='black', alpha=0.7, density=True, histtype="step", label='T = 5')  # Histograma ponderado por Pi_m
# plt.grid(alpha=.6, linestyle='--')
plt.legend()
# plt.savefig("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/latex/p2/b/N4_energy_GrayFlip.pdf")
plt.show()
