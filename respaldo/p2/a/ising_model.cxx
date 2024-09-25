#include "IsingModel.h"
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>

// Constructor de la clase
IsingModel::IsingModel(int N) : N(N) {
    spins = createSpinsGrid();
}

// Función para crear una red de spins aleatorios (+1 o -1)
std::vector<std::vector<int>> IsingModel::createSpinsGrid() {
    std::vector<std::vector<int>> grid(N, std::vector<int>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            grid[i][j] = rand() % 2 == 0 ? 1 : -1;  // Asigna +1 o -1
        }
    }
    return grid;
}

// Función para imprimir la configuración de spins en formato de matriz
void IsingModel::printSpins() const {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << (spins[i][j] == 1 ? "+" : "-") << " ";
        }
        std::cout << std::endl;
    }
}

// Función para calcular la energía con condiciones de borde periódicas
int IsingModel::energyIsingPeriodic(int J, const std::vector<std::vector<int>>& spins) const {
    int energy = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int right = (i + 1) % N;  // Vecino a la derecha
            int down = (j + 1) % N;   // Vecino de abajo

            energy -= J * spins[i][j] * (spins[right][j] + spins[i][down]);
        }
    }
    return energy;
}

// Función para calcular la energía con condiciones de borde no periódicas
int IsingModel::energyIsingNonPeriodic(int J, const std::vector<std::vector<int>>& spins) const {
    int energy = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i + 1 < N) {
                energy -= J * spins[i][j] * spins[i + 1][j];  // Vecino a la derecha
            }
            if (j + 1 < N) {
                energy -= J * spins[i][j] * spins[i][j + 1];  // Vecino de abajo
            }
        }
    }
    return energy;
}

// Función para enumerar todas las configuraciones posibles de spins
std::vector<std::vector<std::vector<int>>> IsingModel::EnumerateIsing() const {
    std::vector<std::vector<std::vector<int>>> configurations;
    long long totalConfigs = 1LL << (N * N);  // Usar long long para evitar desbordamiento

    for (long long config = 0; config < totalConfigs; ++config) {
        std::vector<int> spinsFlat(N * N);
        for (int i = 0; i < N * N; ++i) {
            spinsFlat[i] = ((config >> i) & 1) == 1 ? 1 : -1;  // 1 para spin up, -1 para spin down.
        }

        // Convertir el vector plano en una matriz de NxN.
        std::vector<std::vector<int>> spinsGrid(N, std::vector<int>(N));
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                spinsGrid[i][j] = spinsFlat[i * N + j];
            }
        }

        configurations.push_back(spinsGrid);
    }

    return configurations;
}



// Función gray-flip para cambiar un solo spin entre configuraciones
void IsingModel::gray_flip(const std::vector<int>& previousGray, const std::vector<int>& currentGray) {
    for (int i = 0; i < N * N; ++i) {
        if (previousGray[i] != currentGray[i]) {
            int row = i / N;
            int col = i % N;
            spins[row][col] *= -1;  // Flip del spin
            break;
        }
    }
}
