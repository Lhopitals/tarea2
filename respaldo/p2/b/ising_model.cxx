#include "IsingModel.h"
#include <iostream>
#include <cstdlib>
#include <ctime>

// Constructor de la clase
IsingModel::IsingModel(int N) : N(N) {
    spins = createSpinsGrid();
}

// Función para generar un número aleatorio uniforme
double IsingModel::random_uniform(double min, double max) {
    return min + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX) / (max - min));
}

// Función para crear una red de spins aleatorios (+1 o -1)
std::vector<std::vector<int>> IsingModel::createSpinsGrid() {
    std::vector<std::vector<int>> grid(N, std::vector<int>(N));
    srand(time(0));  // Semilla para aleatoriedad

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
            if (spins[i][j] == -1) {
                std::cout << "- ";  // Representa -1 como "-"
            } else {
                std::cout << "+ ";  // Representa 1 como "+"
            }
        }
        std::cout << std::endl;
    }
}

int IsingModel::energyIsingPeriodic(int N, int J, const std::vector<std::vector<int>>& spins) const {
    int energy = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int right = (i + 1) % N;  // Vecino a la derecha
            int down = (j + 1) % N;   // Vecino de abajo

            // Sumar interacciones de los spins vecinos
            energy -= J * spins[i][j] * (spins[right][j] + spins[i][down]);
        }
    }
    return energy;
}


// Función para calcular la energía con condiciones de borde no periódicas
int IsingModel::energyIsingNonPeriodic(int N, int J, const std::vector<std::vector<int> >& spins) const {
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

// Getter para obtener los spins
const std::vector<std::vector<int>>& IsingModel::getSpins() const {
    return spins;
}

// Setter para actualizar los spins
void IsingModel::setSpins(const std::vector<std::vector<int>>& newSpins) {
    spins = newSpins;
}

// // Función de Gray Flip que cambia un solo spin
// void IsingModel::gray_flip(const std::vector<int>& previousGray, const std::vector<int>& currentGray) {
//     // Detectar el índice del spin que ha cambiado comparando el vector actual con el anterior
//     for (int i = 0; i < N * N; ++i) {
//         if (previousGray[i] != currentGray[i]) {
//             // Convertir el índice lineal 'i' a coordenadas 2D (fila y columna) para la matriz de spins
//             int row = i / N;
//             int col = i % N;

//             // Realizar el flip del spin en la posición correspondiente
//             spins[row][col] *= -1;  // Flip en el índice que cambió
//             break;
//         }
//     }
// }

std::vector<std::vector<std::vector<int>>> IsingModel::EnumerateIsing(int N)  const {
    std::vector<std::vector<std::vector<int>>> configurations;  // Vector de matrices

    int totalConfigs = 1 << (N * N);  // Hay 2^(N*N) configuraciones posibles.

    for (int config = 0; config < totalConfigs; ++config) {
        std::vector<int> spinsFlat(N * N);  // Vector plano de tamaño N*N.

        // Convertimos el número `config` en una configuración de spins.
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

        // Agregar la configuración (matriz de NxN) a la lista de configuraciones.
        configurations.push_back(spinsGrid);
    }

    return configurations;  // Devolver el vector de matrices
}

// Función para calcular la función partición Z(β)
double IsingModel::PartitionFunction(int J, double beta) const {
    // int totalConfigs = 1 << (N * N);  // Número total de configuraciones (2^(N*N))
    double partitionFunction = 0.0;

    // Enumerar todas las configuraciones posibles
    std::vector<std::vector<std::vector<int>>> configs = EnumerateIsing(N);

    // Calcular Z(β)
    for (const auto& spins : configs) {
        int energy = energyIsingPeriodic(N, J, spins);
        partitionFunction += exp(-beta * energy);  // Agregar contribución de cada configuración
    }

    return partitionFunction;
}

//Calculo la energía promedio
double IsingModel::AverageEnergy(int J, double beta, int samplingSteps) {
    int energy = energyIsingPeriodic(N, J, spins);
    double totalEnergy = 0.0;

    // Realizamos la simulación de Metropolis y acumulamos la energía
    for (int step = 0; step < samplingSteps; ++step) {
        metropolis(N, J, 1.0 / beta, energy);  // Actualizamos el sistema usando Metropolis
        totalEnergy += energy;
    }

    // Promediamos la energía
    return totalEnergy / samplingSteps;
}

//Calculo la energía promedio
double IsingModel::SquareAverageEnergy(int J, double beta, int samplingSteps) {
    int energy = energyIsingPeriodic(N, J, spins);
    double totalEnergySquared = 0.0;

    // Realizamos la simulación de Metropolis y acumulamos la energía cuadrada
    for (int step = 0; step < samplingSteps; ++step) {
        metropolis(N, J, 1.0 / beta, energy);  // Actualizamos el sistema usando Metropolis
        totalEnergySquared += energy * energy;
    }

    // Promediamos la energía cuadrada
    return totalEnergySquared / samplingSteps;
}

// Función para calcular la magnetización
double IsingModel::magnetization() const {
    int totalMagnetization = 0;

    // Sumar todos los spins
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            totalMagnetization += spins[i][j];
        }
    }

    // Devolver la magnetización absoluta promedio por sitio
    return std::abs(totalMagnetization) / static_cast<double>(N * N);
}

// Función para generar el código Gray para un número (64 bits)
unsigned long long grayCode(unsigned long long num) {
    return num ^ (num >> 1);  // Regresa el código Gray para el número
}


// Función para hacer flip de un solo spin usando Gray codes
void IsingModel::gray_flip(const std::vector<int>& previousGray, const std::vector<int>& currentGray) {
    for (int i = 0; i < N * N; ++i) {
        if (previousGray[i] != currentGray[i]) {
            int row = i / N;
            int col = i % N;
            spins[row][col] *= -1;  // Flip del spin en la posición que cambió
            break;
        }
    }
}

// Función para contar configuraciones con energías y magnetizaciones específicas
std::map<std::pair<int, int>, int> IsingModel::countGrayConfigurationsEnergyMagnetization(int N, int J) {
    std::map<std::pair<int, int>, int> histogram;  // (E, M) -> count
    unsigned long long totalConfigs = 1ULL << (N * N);  // 2^(N*N) configuraciones posibles

    // Configuración inicial: todos los spins +1
    std::vector<int> spinsFlat(N * N, 1);  // Comienza con todos los spins up
    std::vector<int> prevGrayCode = spinsFlat;
    setSpins(spinsFlatToGrid(spinsFlat));  // Asegura que 'spins' coincide con 'spinsFlat'

    for (unsigned long long config = 0; config < totalConfigs; ++config) {
        // Generar código Gray
        unsigned long long gray = grayCode(config);
        std::vector<int> currentGray(N * N);
        for (int i = 0; i < N * N; ++i) {
            currentGray[i] = (gray & (1ULL << i)) ? 1 : -1;  // Convertimos de Gray a spin
        }

        // Usar Gray flip para cambiar solo un spin entre configuraciones consecutivas
        gray_flip(prevGrayCode, currentGray);
        prevGrayCode = currentGray;

        // Calcular energía y magnetización
        int E = energyIsingPeriodic(N, J, spins);
        int M = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                M += spins[i][j];
            }
        }

        histogram[std::make_pair(E, M)]++;  // Incrementar el conteo para esa combinación de (E, M)
    }

    return histogram;
}


// Función auxiliar para convertir spinsFlat a una matriz 2D
std::vector<std::vector<int>> IsingModel::spinsFlatToGrid(const std::vector<int>& spinsFlat) const {
    std::vector<std::vector<int>> spinsGrid(N, std::vector<int>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            spinsGrid[i][j] = spinsFlat[i * N + j];
        }
    }
    return spinsGrid;
}



std::map<int, double> IsingModel::magnetizationDistribution(double T, int J) {
    std::map<int, double> probM;
    double Z = 0.0;  // Función partición

    std::map<std::pair<int, int>, int> hist = countGrayConfigurationsEnergyMagnetization(N, J);

    for (const auto& entry : hist) {
        int E = entry.first.first;
        int M = entry.first.second;
        int count = entry.second;

        double weight = count * exp(-E / T);
        probM[M] += weight;
        Z += weight;
    }

    for (auto& entry : probM) {
        entry.second /= Z;  // Normalizar para obtener la probabilidad
    }

    return probM;
}


double IsingModel::binderCumulant(double T, int J) {
    std::map<int, double> probM = magnetizationDistribution(T, J);
    double m2 = 0.0, m4 = 0.0;

    for (const auto& entry : probM) {
        double m = entry.first / static_cast<double>(N * N);  // Magnetización por spin
        m2 += entry.second * m * m ;
        m4 += entry.second * m * m * m * m;
    }

    return 0.5 * (3.0 - m4 / (m2 * m2));
}

// Función para calcular el cambio de energía si el spin (i, j) se voltea
int IsingModel::deltaE(int N, int i, int j, int J) {
    int up = (i - 1 + N) % N;
    int down = (i + 1) % N;
    int left = (j - 1 + N) % N;
    int right = (j + 1) % N;

    int sumNeighbors = spins[up][j] + spins[down][j] + spins[i][left] + spins[i][right];
    return 2 * J * spins[i][j] * sumNeighbors;  // Diferencia de energía
}


void IsingModel::metropolis(int N, int J, double T, int &E) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int delta_E = deltaE(N, i, j, J);

            // Si el cambio de energía es negativo o aceptamos el cambio por Metropolis
            if (delta_E <= 0 || random_uniform(0.0, 1.0) < exp(-delta_E / T)) {
                spins[i][j] *= -1;  // Voltea el spin
                E += delta_E;       // Actualiza la energía total del sistema
            }
        }
    }
}



