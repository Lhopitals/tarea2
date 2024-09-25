#include "IsingModel.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <map>

int main() {
    // int N = 6;  // Tamaño de la red
    // int J = 1;   // Constante de interacción
    // double T_histogram = 5.0;  // Temperatura fija para el histograma
    // double T_min = 0.1;  // Temperatura mínima para el Binder cumulant (evitar T=0)
    // double T_max = 5.0;  // Temperatura máxima para el Binder cumulant
    // double dT = 0.1;     // Incremento de temperatura
    // int equilibrationSteps = 1e5;  // Pasos de equilibración
    // int samplingSteps = 1e6;  // Pasos para tomar muestras

    // // Archivos de salida
    // // std::ofstream binderFile("N32_binder_vs_temperatureT5.dat");
    // std::ofstream magnetizationFile("N6_magnetization_histogram_T5.dat");

    // // Crear un objeto de la clase IsingModel
    // IsingModel ising(N);

    // // Fase 1: Calcular el histograma de la magnetización para una temperatura fija (T = T_histogram)
    // int energy = ising.energyPeriodic(J);

    // // Fase de equilibración con Wolff Cluster
    // for (int step = 0; step < equilibrationSteps; ++step) {
    //     ising.wolffCluster(J, T_histogram, energy);
    // }

    // // Variables para almacenar el histograma de magnetización
    // std::map<int, int> magnetization_histogram;

    // // Calcular la magnetización para T = T_histogram
    // for (int step = 0; step < samplingSteps; ++step) {
    //     ising.wolffCluster(J, T_histogram, energy);

    //     // Calcular la magnetización total (incluyendo valores negativos)
    //     int M = static_cast<int>(ising.magnetization() * N * N);
    //     magnetization_histogram[M]++;
    // }

    // // Guardar el histograma de magnetización en un archivo
    // for (const auto& entry : magnetization_histogram) {
    //     magnetizationFile << entry.first << " " << entry.second << std::endl;
    // }

    // magnetizationFile.close();

    // // Fase 2: Calcular el Binder cumulant para una lista de temperaturas entre T_min y T_max
    // for (double T = T_min; T <= T_max; T += dT) {
    //     // Inicializar energía para cada temperatura
    //     energy = ising.energyPeriodic(J);

    //     // Fase de equilibración con Wolff Cluster
    //     for (int step = 0; step < equilibrationSteps; ++step) {
    //         ising.wolffCluster(J, T, energy);
    //     }

    //     // Variables para calcular el Binder cumulant
    //     double avgM2 = 0.0, avgM4 = 0.0;

    //     // Calculo de la magnetización para el Binder cumulant
    //     for (int step = 0; step < samplingSteps; ++step) {
    //         ising.wolffCluster(J, T, energy);

    //         // Acumular para el Binder cumulant
    //         double m = ising.magnetization();
    //         avgM2 += m * m;
    //         avgM4 += m * m * m * m;
    //     }

    //     // Promediar los valores acumulados
    //     avgM2 /= samplingSteps;
    //     avgM4 /= samplingSteps;

    //     // Calcular el Binder cumulant
    //     double binder = (0.5) * (3.0 - (avgM4 / (avgM2 * avgM2)));

    //     // Guardar resultados del Binder cumulant
    //     binderFile << std::fixed << std::setprecision(7) << T << " " << binder << std::endl;

    //     // Mostrar el resultado en la terminal
    //     std::cout << "T: " << T << " Binder: " << binder << std::endl;
    // }

    // binderFile.close();
    
int N = 6;  // Tamaño de la red
int J = 1;   // Constante de interacción
int equilibrationSteps = static_cast<int>(2e5);  // Pasos de equilibración
int samplingSteps = static_cast<int>(5e6);  // Pasos de muestreo
double T_fija = 2.0;  // Temperatura fija para análisis detallado
double beta_fija = 1.0 / T_fija;

// Crear un objeto de la clase IsingModel para red de tamaño N
IsingModel ising(N);
std::cout << "Red de tamaño N = " << N << " creada\n";

// **Primera Parte: Cálculo de Energía y Calor Específico para T fija**
int energy = ising.energyPeriodic(J);  // Cálculo inicial de la energía usando condiciones periódicas

// Fase de equilibración con el algoritmo Wolff
for (int step = 0; step < equilibrationSteps; ++step) {
    ising.wolffCluster(J, T_fija, energy);  // Actualizamos la configuración usando el método Wolff
}

for (int i = 0; i < 5; ++i) {
    // Calcular energía y calor específico usando el algoritmo Wolff
double E_mean = ising.AverageEnergyWolff(J, beta_fija, samplingSteps);  // Energía promedio
double E2_mean = ising.SquareAverageEnergyWolff(J, beta_fija, samplingSteps);  // Energía cuadrada promedio
double heat_capacity = (E2_mean - E_mean * E_mean) / (T_fija * T_fija * N * N);  // Calor específico

std::cout << "T: " << T_fija << ", <E/N>: " << E_mean / (N * N)
          << ", C_v: " << heat_capacity << std::endl;

}
    return 0;
}
