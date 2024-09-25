#include "IsingModel.h"
#include <iostream>
#include <vector>
#include <map>
#include <cmath>  // Para exp()
#include <fstream>  // Para escribir archivos

int main() {
    int N = 2;  // Tamaño de la red (puedes cambiarlo a 2, 4 o 6)
    int J = 1;  // Constante de interacción
    double T = 2.0;  // Temperatura
    double beta = 1.0 / T;  // Beta = 1/T

    // Crear un objeto de la clase IsingModel
    IsingModel ising(N);

    // Mapas para contar configuraciones con energías y magnetizaciones específicas
    std::map<std::pair<int, int>, int> energyMagnetizationCounts;  // (E, M) -> número de configuraciones
    std::map<int, int> magnetizationHistogram;  // M -> número de configuraciones
    std::map<double, double> magnetizationDistribution;  // m -> distribución de probabilidad

    // Archivos de salida
    std::ofstream file_histograma_EM("N2_N(E,M)T2.dat");
    std::ofstream file_histograma_E("N2_E_GrayFlipT2.dat");
    std::ofstream file_histograma_M("N2_M_GrayFlipT2.dat");  // Histograma de magnetización
    std::ofstream file_distribucion_M("N2_magnetizacion_distribucion_GrayFlipT2.dat");
    std::ofstream file_binder("N2_binder_cumulant_GrayFlip.dat");

    // Usar el algoritmo Gray Flip para contar configuraciones con energía y magnetización
    energyMagnetizationCounts = ising.countGrayConfigurationsEnergyMagnetization(N, J);

    // Escribir el histograma N(E, M) a un archivo .dat
    // file_histograma_EM << "# Energía\tMagnetización\tCantidad\n";
    for (const auto& entry : energyMagnetizationCounts) {
        file_histograma_EM << entry.first.first << "\t" << entry.first.second << "\t" << entry.second << "\n";
    }
    file_histograma_EM.close();
    std::cout << "Histograma N(E, M) escrito usando Gray Flip\n";

    // Escribir el histograma N(E) (sumado sobre M) a un archivo .dat con suma por energía
    // file_histograma_E << "Energía\tCantidad Sumada\n";
    std::map<int, int> energySum;

    // Sumar las configuraciones con la misma energía
    for (const auto& entry : energyMagnetizationCounts) {
        int energy = entry.first.first;
        energySum[energy] += entry.second;
    }

    // Escribir los resultados acumulados en el archivo
    for (const auto& entry : energySum) {
        file_histograma_E << entry.first << "\t" << entry.second << "\n";
    }
    file_histograma_E.close();
    std::cout << "Histograma N(E) escrito usando Gray Flip\n";

    // Calcular el histograma de la magnetización
    for (const auto& entry : energyMagnetizationCounts) {
        int magnetization = entry.first.second;
        magnetizationHistogram[magnetization] += entry.second;
    }

    // Escribir el histograma de magnetización N(M) a un archivo .dat
    // file_histograma_M << "# Magnetización\tCantidad\n";
    for (const auto& entry : magnetizationHistogram) {
        file_histograma_M << entry.first << "\t" << entry.second << "\n";
    }
    file_histograma_M.close();
    std::cout << "Histograma N(M) escrito usando Gray Flip\n";

    // Calcular la distribución de probabilidad de la magnetización por sitio, m = M / N^2
    // file_distribucion_M << "# Magnetización por sitio m\tProbabilidad\n";
    for (int M = -N * N; M <= N * N; M += 2) {
        double m = static_cast<double>(M) / (N * N);  // Magnetización por sitio
        double prob = 1.0 / (1.0 + std::exp(-2.0 * beta * m));
        magnetizationDistribution[m] = prob;
        file_distribucion_M << m << "\t" << prob << "\n";
    }
    file_distribucion_M.close();
    std::cout << "Distribución de magnetización por sitio π(m) escrita\n";

    double T_min = 0.1;  // Temperatura mínima
    double T_max = 5;  // Temperatura máxima
    double dT = 0.1;  // Incremento de temperatura

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
