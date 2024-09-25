import numpy as np
import matplotlib.pyplot as plt

# Uso de Latex en los ejes
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams['font.size'] = 14
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "serif"

# Leer archivo T = 2
with open("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p4/N6_magnetization_histogram_T2.dat", "r") as f1:
    lines1 = f1.readlines()
    M1 = [float(line.split()[0]) for line in lines1]  # Magnetización (primera columna)
    count1 = [float(line.split()[1]) for line in lines1]  # Probabilidad (segunda columna)

# Leer archivo para T = 5
with open("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p4/N6_magnetization_histogram_T5.dat", "r") as f2:
    lines2 = f2.readlines()
    M2 = [float(line.split()[0]) for line in lines2]  # Magnetización (primera columna)
    count2 = [float(line.split()[1]) for line in lines2]  # Probabilidad (segunda columna)

# Agregar los valores de magnetización negativos, excepto para el 0
M1_extended = M1 + [-m for m in M1 if m != 0]
count1_extended = count1 + [c for m, c in zip(M1, count1) if m != 0]  # Duplicar solo para los valores no nulos

M2_extended = M2 + [-m for m in M2 if m != 0]
count2_extended = count2 + [c for m, c in zip(M2, count2) if m != 0]  # Duplicar solo para los valores no nulos

# Crear histograma ponderado por las probabilidades
plt.xlabel("$M$")
plt.hist(M1_extended, bins=21, weights=count1_extended, edgecolor='black', alpha=0.7, density=True, histtype="step", label='T = 2', linestyle="--")  # Histograma ponderado por Pi_m
plt.hist(M2_extended, bins=21, weights=count2_extended, edgecolor='black', alpha=0.7, density=True, histtype="step", label='T = 5')  # Histograma ponderado por Pi_m
# plt.savefig
plt.legend()
plt.show()
