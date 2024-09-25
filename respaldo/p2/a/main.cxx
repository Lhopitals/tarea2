#include "IsingModel.h"
#include <iostream>
#include <fstream>  // Para escribir en archivos
#include <vector>
#include <map>

int main() {
    int N = 6;  // Cambia el tamaño de la red a 6x6
    int J = 1;  // Constante de interacción
    long long blockSize = 1LL << 30;  // Procesar en bloques de configuraciones

    // Crear un objeto de la clase IsingModel
    IsingModel ising(N);

    std::map<int, long long> energyCountsNonPeriodic;  // Usar long long para grandes números (no periódicas)
    std::map<int, long long> energyCountsPeriodic;     // Usar long long para grandes números (periódicas)

    // Abrir archivos para escribir
    std::ofstream nonPeriodicFile("ising_non_periodic_results.dat");
    std::ofstream periodicFile("ising_periodic_results.dat");

    if (!nonPeriodicFile.is_open() || !periodicFile.is_open()) {
        std::cerr << "Error al abrir los archivos para escribir los resultados." << std::endl;
        return 1;  // Salir si no se pueden abrir los archivos
    }

    for (long long blockStart = 0; blockStart < (1LL << (N * N)); blockStart += blockSize) {
        long long blockEnd = std::min(blockStart + blockSize, 1LL << (N * N));

        for (long long config = blockStart; config < blockEnd; ++config) {
            // Crear la configuración usando bits
            std::vector<std::vector<int>> spinsGrid(N, std::vector<int>(N));
            for (int i = 0; i < N * N; ++i) {
                spinsGrid[i / N][i % N] = ((config >> i) & 1) ? 1 : -1;
            }

            // Calcular la energía de la configuración con condiciones no periódicas
            int energyNonPeriodic = ising.energyIsingNonPeriodic(J, spinsGrid);
            energyCountsNonPeriodic[energyNonPeriodic]++;

            // Calcular la energía de la configuración con condiciones periódicas
            int energyPeriodic = ising.energyIsingPeriodic(J, spinsGrid);
            energyCountsPeriodic[energyPeriodic]++;
        }
    }

    // Escribir los resultados de condiciones de borde no periódicas en el archivo
    nonPeriodicFile << "Condiciones de borde: No Periódicas." << std::endl;
    nonPeriodicFile << "Tamaño de la red: " << N << "x" << N << std::endl;
    for (const auto& pair : energyCountsNonPeriodic) {
        nonPeriodicFile << "Cantidad de configuraciones con energía " << pair.first << ": " << pair.second << std::endl;
    }

    // Escribir los resultados de condiciones de borde periódicas en el archivo
    periodicFile << "Condiciones de borde: Periódicas." << std::endl;
    periodicFile << "Tamaño de la red: " << N << "x" << N << std::endl;
    for (const auto& pair : energyCountsPeriodic) {
        periodicFile << "Cantidad de configuraciones con energía " << pair.first << ": " << pair.second << std::endl;
    }

    // Cerrar los archivos
    nonPeriodicFile.close();
    periodicFile.close();

    std::cout << "Los resultados han sido guardados en 'ising_non_periodic_results.dat' y 'ising_periodic_results.dat'." << std::endl;

    return 0;
}
