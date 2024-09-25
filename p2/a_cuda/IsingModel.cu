#include "IsingModel.h"
#include <iostream>
#include <cuda.h>
#include <cmath>
#include <vector>

// Kernel CUDA para calcular las energías
__global__ void calculateEnergiesKernel(int* spins, int N, int J, int* energiesNonPeriodic, int* energiesPeriodic, int totalConfigs) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= totalConfigs) return;

    int row, col, right, down;
    int energyNonPeriodic = 0, energyPeriodic = 0;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            row = i;
            col = j;
            right = (i + 1) % N;  // Vecino a la derecha con condición periódica
            down = (j + 1) % N;   // Vecino de abajo con condición periódica

            // Energía no periódica
            if (i + 1 < N) {
                energyNonPeriodic -= J * spins[idx * N * N + i * N + j] * spins[idx * N * N + (i + 1) * N + j];
            }
            if (j + 1 < N) {
                energyNonPeriodic -= J * spins[idx * N * N + i * N + j] * spins[idx * N * N + i * N + (j + 1)];
            }

            // Energía periódica
            energyPeriodic -= J * spins[idx * N * N + i * N + j] * (spins[idx * N * N + right * N + j] + spins[idx * N * N + i * N + down]);
        }
    }

    // Guardar resultados en arrays globales
    energiesNonPeriodic[idx] = energyNonPeriodic;
    energiesPeriodic[idx] = energyPeriodic;
}

// Constructor
IsingModel::IsingModel(int N) : N(N) {
    // Inicializa las configuraciones de spins usando EnumerateIsing
    spins = createSpinsGrid();
}

// Destructor
IsingModel::~IsingModel() {
    // Nada que liberar en este caso
}

// Método para enumerar todas las configuraciones posibles de spins
std::vector<std::vector<std::vector<int>>> IsingModel::EnumerateIsing() const {
    std::vector<std::vector<std::vector<int>>> configurations;
    long long totalConfigs = 1LL << (N * N);  // 2^(N*N) configuraciones

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

// Función que calcula las energías usando enumeración de Ising y CUDA
void IsingModel::calculateEnergiesWithEnumerateIsingCUDA(int J, std::vector<int>& energiesNonPeriodic, std::vector<int>& energiesPeriodic) {
    std::vector<std::vector<std::vector<int>>> configurations = EnumerateIsing();
    long long totalConfigs = configurations.size();

    // Alocar memoria para almacenar spins en formato plano
    int* spinsFlat = new int[totalConfigs * N * N];
    for (long long config = 0; config < totalConfigs; ++config) {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                spinsFlat[config * N * N + i * N + j] = configurations[config][i][j];
            }
        }
    }

    // Alocar memoria en GPU
    int *d_spins, *d_energiesNonPeriodic, *d_energiesPeriodic;
    cudaMalloc(&d_spins, totalConfigs * N * N * sizeof(int));
    cudaMalloc(&d_energiesNonPeriodic, totalConfigs * sizeof(int));
    cudaMalloc(&d_energiesPeriodic, totalConfigs * sizeof(int));

    // Copiar spins a la GPU
    cudaMemcpy(d_spins, spinsFlat, totalConfigs * N * N * sizeof(int), cudaMemcpyHostToDevice);

    // Configurar tamaño de bloques e hilos
    int blockSize = 256;
    int numBlocks = (totalConfigs + blockSize - 1) / blockSize;

    // Lanzar kernel para calcular energías
    calculateEnergiesKernel<<<numBlocks, blockSize>>>(d_spins, N, J, d_energiesNonPeriodic, d_energiesPeriodic, totalConfigs);

    // Copiar los resultados a la CPU
    cudaMemcpy(energiesNonPeriodic.data(), d_energiesNonPeriodic, totalConfigs * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(energiesPeriodic.data(), d_energiesPeriodic, totalConfigs * sizeof(int), cudaMemcpyDeviceToHost);

    // Liberar memoria
    cudaFree(d_spins);
    cudaFree(d_energiesNonPeriodic);
    cudaFree(d_energiesPeriodic);
    delete[] spinsFlat;
}
