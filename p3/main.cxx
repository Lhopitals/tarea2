#include "IsingModel.h"
#include <iostream>
#include <vector>
#include <fstream> 
#include <cmath>

int main() {
    int N = 6;  // Tamaño de la red
    int J = 1;  // Constante de interacción
    double T = 2.0;  // Temperatura
    int equilibrationSteps = 1e5;  // Pasos de equilibración
    int SamplingSteps = 10e6;


    // Crear un objeto de la clase IsingModel
    IsingModel ising(N);

    // Fase de equilibración: No tomo
    int energy = ising.energyIsingPeriodic(N, J, ising.getSpins());
    for (int step = 0; step < equilibrationSteps; ++step) {
        ising.metropolis(N, J, T, energy);
    }

    double beta = 1.0 / T;  // Inverso de la temperatura
    for (int i = 0; i < 5; i++){
            // Utilizar AverageEnergy y SquareAverageEnergy
            double E_mean = ising.AverageEnergy(J, beta, SamplingSteps);           // Energía promedio
            double E2_mean = ising.SquareAverageEnergy(J, beta, SamplingSteps);    // Energía cuadrática promedio
            double heat_capacity = (E2_mean - E_mean * E_mean) / (T * T * N * N);

            // Imprimir resultados
            std::cout << "Temperatura : " << T << ", " <<" < E/N > : " << E_mean / (N * N) << ", " << " C_v : " << heat_capacity << std::endl;
        }
    //      // Calcula la magnetización promedio
    // double mag = ising.magnetization();
    // std::cout << "Magnetización promedio por sitio: " << mag << std::endl;

    
    // int N2 = 32;  // Tamaño de la red
    // IsingModel ising2(N2);
    // double T_min = 1.0;  // Temperatura mínima
    // double T_max = 5.0;  // Temperatura máxima
    // int num_T = 100;  // Número de puntos de temperatura
    // // int equilibrationSteps = 1e5;  // Pasos de equilibración
    // // int samplingSteps = 1e6;  // Pasos para tomar muestras

    // // Vector para almacenar la magnetización promedio por temperatura
    // std::vector<double> magnetizations;

    // // Archivo para guardar los resultados
    // std::ofstream outfile("N32_magnetization_vs_temperature.dat");

    // // Bucle sobre diferentes temperaturas
    // for (int t = 0; t < num_T; ++t) {
    //     double T = T_min + t * (T_max - T_min) / (num_T - 1);
    //     double beta = 1.0 / T;

    //     // Inicializar energía
    //     int energy2 = ising2.energyIsingPeriodic(N2, J, ising2.getSpins());

    //     // Fase de equilibración
    //     for (int step = 0; step < equilibrationSteps; ++step) {
    //         ising2.metropolis(N2, J, T, energy);
    //     }

    //     // Medición de la magnetización
    //     double magnetization_sum = 0.0;

    //     // Realizar mediciones
    //     for (int step = 0; step < samplingSteps; ++step) {
    //         ising2.metropolis(N2, J, T, energy2);
    //         magnetization_sum += ising2.magnetization();
    //     }

    //     // Promediar la magnetización
    //     double magnetization_mean = magnetization_sum / samplingSteps;
    //     magnetizations.push_back(magnetization_mean);

    //     // Guardar los resultados en el archivo
    //     outfile << T << " " << magnetization_mean << std::endl;

    //     // Imprimir resultados
    //     std::cout << "T: " << T << " Magnetización promedio por sitio: " << magnetization_mean << std::endl;
    // }

    // outfile.close();  // Cerrar el archivo

    
    return 0;
}
