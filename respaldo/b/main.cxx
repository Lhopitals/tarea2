#include "IsingModel.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>  // Para exp()
#include <fstream>  // Para escribir archivos

int main() {
    int N = 2;  // Tamaño de la red (puedes cambiarlo a 2, 4 o 6)
    int J = 1;  // Constante de interacción
    double T = 2.5;  // Temperatura (puedes variar la temperatura para estudiar el sistema)
    double beta = 1.0 / T;  // Beta = 1/T

    // Crear un objeto de la clase IsingModel
    IsingModel ising(N);

    // Mapas para contar configuraciones con energías y magnetizaciones específicas
    std::map<std::pair<int, int>, int> energyMagnetizationCounts;  // (E, M) -> número de configuraciones
    std::map<int, int> energyCounts;  // E -> número total de configuraciones (sumando sobre M)

    // Archivos de salida
    std::ofstream file_histograma_EM("N2_N(E,M)_GrayFlipT2.dat");
    std::ofstream file_histograma_E("N2_E_GrayFlipT2.dat");
    std::ofstream file_distribucion_M("N2_magnetizacion_GrayFlipT2.dat");
    std::ofstream file_binder("N2_binder_cumulant_GrayFlipT2.dat");

    // Usar el algoritmo Gray Flip para contar configuraciones con energía y magnetización
    energyMagnetizationCounts = ising.countGrayConfigurationsEnergyMagnetization(N, J);

    // Escribir el histograma N(E, M) a un archivo .dat
    file_histograma_EM << "# Energía\tMagnetización\tCantidad\n";
    for (const auto& entry : energyMagnetizationCounts) {
        file_histograma_EM << entry.first.first << "\t" << entry.first.second << "\t" << entry.second << "\n";
    }
    file_histograma_EM.close();
    std::cout << "Histograma N(E, M) escrito usando Gray Flip\n";

    // Escribir el histograma N(E) (sumado sobre M) a un archivo .dat con suma por energía
    file_histograma_E << "Energía\tCantidad Sumada\n";
    std::map<int, int> energySum;  // Mapa para almacenar la suma de las configuraciones por energía

    // Sumar las configuraciones con la misma energía
    for (const auto& entry : energyMagnetizationCounts) {
        int energy = entry.first.first;  // Extraer la energía
        energySum[energy] += entry.second;
        // if (energy >= 0) {  // Dejar pasar solo las energías positivas
        //     energySum[energy] += entry.second;  // Sumar la cantidad para esa energía
        // }
    }

    // Escribir los resultados acumulados en el archivo
    for (const auto& entry : energySum) {
        file_histograma_E << entry.first << "\t" << entry.second << "\n";  // Energía y suma acumulada
    }

    file_histograma_E.close();

    std::cout << "Histograma N(E) escrito usando Gray Flip\n";

    // Calcular la distribución de probabilidad de la magnetización, π_M
    std::map<int, double> magnetizationDistribution;
    double Z = 0.0;  // Función partición

    for (const auto& entry : energyMagnetizationCounts) {
        int energy = entry.first.first;
        int magnetization = entry.first.second;

        double weight = entry.second * exp(-beta * energy);  // e^(-beta * E)
        magnetizationDistribution[magnetization] += weight;
        Z += weight;
    }

    // Normalizar la distribución de magnetización
    for (auto& entry : magnetizationDistribution) {
        entry.second /= Z;  // Dividir por la función partición para normalizar
    }

    // Escribir la distribución de magnetización π(M) a un archivo .dat
    for (const auto& entry : magnetizationDistribution) {
        file_distribucion_M << entry.first << "\t" << entry.second << "\n";
    }
    file_distribucion_M.close();
    std::cout << "Distribución de magnetización π(M) escrita usando Gray Flip \n";

    double T_min = 0.5;  // Temperatura mínima
    double T_max = 100;  // Temperatura máxima
    double dT = 0.5;  // Incremento de temperatura

    // Loop sobre diferentes temperaturas
    for (double T = T_min; T <= T_max; T += dT) {
        // Calcular el Binder cumulant 
        double binderCumulant = ising.binderCumulant(T, J);

        // Escribir el Binder cumulant al archivo
        file_binder << T << "\t" << binderCumulant << "\n";
        std::cout << "T = " << T << ", Binder Cumulant = " << binderCumulant << std::endl;
    }

    file_binder.close();
    std::cout << "Binder Cumulant B(T) escrito usando Gray Flip.\n";

    return 0;
}
