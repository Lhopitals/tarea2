#include "IsingModel.h"

// Constructor de la clase
IsingModel::IsingModel(int N) : N(N) {
    initializeSpins();
}

IsingModel::~IsingModel() {}

// Inicializa la red de spins aleatoriamente
void IsingModel::initializeSpins() {
    spins = std::vector<std::vector<int>>(N, std::vector<int>(N));
    srand(static_cast<unsigned int>(time(nullptr)));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            spins[i][j] = rand() % 2 == 0 ? 1 : -1;
        }
    }
}

// Imprime la configuración de los spins
void IsingModel::printSpins() const {
    for (const auto& row : spins) {
        for (const auto& spin : row) {
            std::cout << (spin == 1 ? "+" : "-") << " ";
        }
        std::cout << std::endl;
    }
}

// Genera un número aleatorio uniforme entre min y max
double IsingModel::randomUniform(double min, double max) {
    return min + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX) / (max - min));
}

// Calcula la energía con condiciones periódicas
int IsingModel::energyPeriodic(int J) const {
    int energy = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int right = (i + 1) % N;
            int down = (j + 1) % N;
            energy -= J * spins[i][j] * (spins[right][j] + spins[i][down]);
        }
    }
    return energy;
}

// Calcula la energía con condiciones no periódicas
int IsingModel::energyNonPeriodic(int J) const {
    int energy = 0;
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i + 1 < N) energy -= J * spins[i][j] * spins[i + 1][j];
            if (j + 1 < N) energy -= J * spins[i][j] * spins[i][j + 1];
        }
    }
    return energy;
}

// Calcula la magnetización total
double IsingModel::magnetization() const {
    int totalMagnetization = 0;
    for (const auto& row : spins) {
        for (int spin : row) {
            totalMagnetization += spin;
        }
    }
    return static_cast<double>(std::abs(totalMagnetization)) / (N * N);
}

// Implementa el algoritmo Metropolis
void IsingModel::metropolis(int J, double T, int& E) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            int delta_E = deltaE(i, j, J);
            if (delta_E <= 0 || randomUniform(0.0, 1.0) < exp(-delta_E / T)) {
                spins[i][j] *= -1;
                E += delta_E;
            }
        }
    }
}

// Implementa el algoritmo Wolff Cluster
void IsingModel::wolffCluster(int J, double T, int& E) {
    double prob = 1 - exp(-2.0 * J / T);
    std::vector<std::vector<bool>> visited(N, std::vector<bool>(N, false));
    int i = rand() % N, j = rand() % N;

    std::vector<std::pair<int, int>> cluster = {{i, j}};
    visited[i][j] = true;
    int spinVal = spins[i][j];

    for (size_t idx = 0; idx < cluster.size(); ++idx) {
        int x = cluster[idx].first;
        int y = cluster[idx].second;
        addToCluster(cluster, visited, x, y, prob, spinVal);
    }

    for (const auto& site : cluster) {
        spins[site.first][site.second] *= -1;
    }
    E = energyPeriodic(J);
}

// Cálculo del cambio de energía al voltear un spin
int IsingModel::deltaE(int i, int j, int J) {
    int up = (i - 1 + N) % N, down = (i + 1) % N;
    int left = (j - 1 + N) % N, right = (j + 1) % N;
    int sumNeighbors = spins[up][j] + spins[down][j] + spins[i][left] + spins[i][right];
    return 2 * J * spins[i][j] * sumNeighbors;
}

// Función privada para obtener los vecinos de un spin
std::vector<std::pair<int, int>> IsingModel::getNeighbors(int i, int j) const {
    return {{(i + 1) % N, j}, {(i - 1 + N) % N, j}, {i, (j + 1) % N}, {i, (j - 1 + N) % N}};
}

// Función para agregar spins vecinos al cluster
void IsingModel::addToCluster(std::vector<std::pair<int, int>>& cluster, std::vector<std::vector<bool>>& visited, int x, int y, double prob, int spinVal) {
    for (const auto& neighbor : getNeighbors(x, y)) {
        int nx = neighbor.first, ny = neighbor.second;
        if (spins[nx][ny] == spinVal && !visited[nx][ny] && randomUniform(0.0, 1.0) < prob) {
            cluster.push_back({nx, ny});
            visited[nx][ny] = true;
        }
    }
}

// Cálculo de la función de partición
double IsingModel::partitionFunction(int J, double beta) const {
    double Z = 0.0;
    int totalConfigs = 1 << (N * N);

    for (int config = 0; config < totalConfigs; ++config) {
        std::vector<int> spinsFlat(N * N);
        for (int i = 0; i < N * N; ++i) {
            spinsFlat[i] = ((config >> i) & 1) == 1 ? 1 : -1;
        }
        int energy = energyPeriodic(J);
        Z += exp(-beta * energy);
    }
    return Z;
}

// Cálculo de la energía promedio usando Metropolis
double IsingModel::averageEnergyMetropolis(int J, double beta, int steps) {
    int energy = energyPeriodic(J);
    double totalEnergy = 0.0;
    for (int step = 0; step < steps; ++step) {
        metropolis(J, 1.0 / beta, energy);
        totalEnergy += energy;
    }
    return totalEnergy / steps;
}

// Cálculo de la energía promedio usando Wolff
double IsingModel::averageEnergyWolff(int J, double beta, int steps) {
    int energy = energyPeriodic(J);
    double totalEnergy = 0.0;
    for (int step = 0; step < steps; ++step) {
        wolffCluster(J, 1.0 / beta, energy);
        totalEnergy += energy;
    }
    return totalEnergy / steps;
}

// Distribución de la magnetización
std::map<int, double> IsingModel::magnetizationDistribution(double T, int J, int steps) {
    std::map<int, int> magCount;
    for (int step = 0; step < steps; ++step) {
        int E = energyPeriodic(J);
        wolffCluster(J, T, E);
        int M = static_cast<int>(magnetization() * N * N);
        magCount[M]++;
    }

    std::map<int, double> distribution;
    for (const auto& [mag, count] : magCount) {
        distribution[mag] = static_cast<double>(count) / steps;
    }
    return distribution;
}

// Cálculo del cumulante de Binder
double IsingModel::binderCumulantMetropolis(double T, int J, int steps) {
    double avgM2 = 0.0, avgM4 = 0.0;
    for (int step = 0; step < steps; ++step) {
        int E = energyPeriodic(J);
        metropolis(J, T, E);
        double M = magnetization();
        avgM2 += M * M;
        avgM4 += M * M * M * M;
    }
    avgM2 /= steps;
    avgM4 /= steps;
    return (0.5)*(3 - (avgM4 / (avgM2 * avgM2)));
}

// Cálculo del cumulante de Binder
double IsingModel::binderCumulantWolff(double T, int J, int steps) {
    double avgM2 = 0.0, avgM4 = 0.0;
    for (int step = 0; step < steps; ++step) {
        int E = energyPeriodic(J);
        wolffCluster(J, T, E);
        double M = magnetization();
        avgM2 += M * M;
        avgM4 += M * M * M * M;
    }
    avgM2 /= steps;
    avgM4 /= steps;
    return (0.5)*(3 - (avgM4 / (avgM2 * avgM2)));
}


// Cálculo de la energía promedio usando el algoritmo Wolff
double IsingModel::AverageEnergyWolff(int J, double beta, int samplingSteps) {
    // Inicializamos la energía con condiciones periódicas
    int energy = energyPeriodic(J);
    double totalEnergy = 0.0;

    // Realizamos la simulación de Wolff y acumulamos la energía
    for (int step = 0; step < samplingSteps; ++step) {
        // Actualizamos el sistema usando el método Wolff
        wolffCluster(J, 1.0 / beta, energy);
        totalEnergy += energy;
    }

    // Promediamos la energía acumulada
    return totalEnergy / samplingSteps;
}

// Cálculo del promedio de la energía cuadrada usando el algoritmo Wolff
double IsingModel::SquareAverageEnergyWolff(int J, double beta, int samplingSteps) {
    // Inicializamos la energía con condiciones periódicas
    int energy = energyPeriodic(J);
    double totalEnergySquared = 0.0;

    // Realizamos la simulación de Wolff y acumulamos la energía cuadrada
    for (int step = 0; step < samplingSteps; ++step) {
        // Actualizamos el sistema usando el método Wolff
        wolffCluster(J, 1.0 / beta, energy);
        totalEnergySquared += energy * energy;
    }

    // Promediamos la energía cuadrada acumulada
    return totalEnergySquared / samplingSteps;
}

