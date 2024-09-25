 #include "IsingModel.h"
#include <iostream>
#include <fstream>
#include <vector>

int main() {
    int N = 6;  // Tamaño de la red
    int J = 1;  // Constante de interacción
    long long totalConfigs = 1LL << (N * N);  // Número total de configuraciones

    // Crear una instancia de la clase IsingModel
    IsingModel ising(N);

    // Alocar memoria para almacenar las energías en la CPU
    std::vector<int> energiesNonPeriodic(totalConfigs);
    std::vector<int> energiesPeriodic(totalConfigs);

    // Calcular energías usando CUDA y EnumerateIsing
    ising.calculateEnergiesWithEnumerateIsingCUDA(J, energiesNonPeriodic, energiesPeriodic);

    // Guardar resultados en archivos
    std::ofstream nonPeriodicFile("ising_non_periodic_cuda_results.dat");
    std::ofstream periodicFile("ising_periodic_cuda_results.dat");

    // Escribir resultados en los archivos
    for (long long i = 0; i < totalConfigs; ++i) {
        nonPeriodicFile << "Configuración " << i << " - Energía No Periódica: " << energiesNonPeriodic[i] << std::endl;
        periodicFile << "Configuración " << i << " - Energía Periódica: " << energiesPeriodic[i] << std::endl;
    }

    // Cerrar los archivos
    nonPeriodicFile.close();
    periodicFile.close();

    std::cout << "Resultados guardados en 'ising_non_periodic_cuda_results.dat' y 'ising_periodic_cuda_results.dat'." << std::endl;

    return 0;
}
