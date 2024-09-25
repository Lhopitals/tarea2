#ifndef ISINGMODEL_H
#define ISINGMODEL_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <map>

class IsingModel {
private:
    int N; // Tamaño de la red
    // int J; // Constante de interacción
    std::vector<std::vector<int>> spins;

public:
    // Constructor de la clase
    IsingModel(int N);

    // Función para generar un número aleatorio uniforme
    double random_uniform(double min, double max);

    // Función para crear una red de spins aleatorios (+1 o -1)
    std::vector<std::vector<int>> createSpinsGrid();

    // Función para imprimir la configuración de spins en formato de matriz
    void printSpins() const;

    // Función para calcular la energía con condiciones de borde periódicas
    int energyIsingPeriodic(int N, int J, const std::vector<std::vector<int> >& spins) const;

    // Función para calcular la energía con condiciones de borde no periódicas
    int energyIsingNonPeriodic(int N, int J, const std::vector<std::vector<int> >& spins) const;

    // Getter para obtener los spins
    const std::vector<std::vector<int>>& getSpins() const;

    // Setter para actualizar los spins
    void setSpins(const std::vector<std::vector<int>>& newSpins);

    //Funcion gray-flip
    void gray_flip(const std::vector<int>& previousGray, const std::vector<int>& currentGray) ;

    //Función Enumerate-Ising
    std::vector<std::vector<std::vector<int>>> EnumerateIsing(int N) const;

    // Función para calcular la función partición Z(β)
    double PartitionFunction(int J, double beta) const;

    // Función para calcular la energía promedio ⟨E⟩
    double AverageEnergy(int J, double beta, int samplingSteps);

    // Función para calcular la energía promedio ⟨E⟩
    double SquareAverageEnergy(int J, double beta, int samplingSteps); 

    // Función para calcular la magnetización promedio
    double magnetization() const;

    std::vector<std::vector<int>> spinsFlatToGrid(const std::vector<int>& spinsFlat) const;

    // Función para calcular la magnetización promedio
    std::map<std::pair<int, int>, int> countGrayConfigurationsEnergyMagnetization(int N, int J);

    // Función para calcular la distribución de probabilidad de la magnetización
    std::map<int, double> magnetizationDistribution(double T, int J);    

    // Función para calcular el Binder Cumulant
    double binderCumulant(double T, int J);


    // Función para calcular el cambio de energía si el spin (i, j) se voltea
    int deltaE(int N, int i, int j, int J);

    // Implementa el algoritmo de Metropolis para una red de spins de tamaño NxN
    void metropolis(int N, int J, double T, int &E);
};

#endif