#!/usr/bin/python3
import math as m
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Uso de LaTeX en los ejes
plt.rcParams.update(plt.rcParamsDefault)
plt.rcParams['font.size'] = 14
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = "serif"

# Leer archivo T = 2
with open("/Users/gabbo/posgrado/magister/cursos/estadistica_avanzada/tarea2/p1/b/rejection_markov.dat", "r") as f1:
    lines1 = f1.readlines()
    N = [float(line.split()[0]) for line in lines1]  # Número de puntos (primera columna)
    Pi_q = [float(line.split()[1]) for line in lines1]  # Error cuadrático medio (segunda columna)

# Configuración de etiquetas de los ejes
plt.xlabel(r"$\delta$")
plt.ylabel(r"$\langle N_{hits}/N - \pi/4\rangle$")
# plt.axhline(y=m.pi/4, color='r', linestyle='--', label=r'$\pi/4$')
plt.plot(N, Pi_q, 'o-', label= "mse-markov-pi")
plt.grid(alpha=.6, linestyle='--')
plt.legend()
plt.show()
