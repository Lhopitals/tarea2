#ifndef ISINGMODEL_H
#define ISINGMODEL_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <map>

class IsingModel {
private:
    int N; // Tamaño de la red
    std::vector<std::vector<int>> spins;

    // Métodos privados para el cálculo de energía
    int deltaE(int i, int j, int J);
    std::vector<std::pair<int, int>> getNeighbors(int i, int j) const;
    void addToCluster(std::vector<std::pair<int, int>>& cluster, std::vector<std::vector<bool>>& visited, int x, int y, double prob, int spinVal);

public:
    // Constructor
    IsingModel(int N);
    ~IsingModel();

    // Métodos públicos
    void initializeSpins();
    void printSpins() const;
    double randomUniform(double min, double max);
    int energyPeriodic(int J) const;
    int energyNonPeriodic(int J) const;
    double magnetization() const;
    void metropolis(int J, double T, int& E);
    void wolffCluster(int J, double T, int& E);

    // Funciones para calcular promedios
    double partitionFunction(int J, double beta) const;
    double averageEnergyMetropolis(int J, double beta, int steps);
    double squareAverageEnergyMetropolis(int J, double beta, int steps);
    double averageEnergyWolff(int J, double beta, int steps);
    double squareAverageEnergyWolff(int J, double beta, int steps);

    // Funciones de análisis
    std::map<std::pair<int, int>, int> countConfigurationsEnergyMagnetization(int J, int steps, double T);
    std::map<int, double> magnetizationDistribution(double T, int J, int steps);
    double binderCumulantMetropolis(double T, int J, int steps);
    double binderCumulantWolff(double T, int J, int steps);
    double AverageEnergyWolff(int J, double beta, int samplingSteps);
    double SquareAverageEnergyWolff(int J, double beta, int samplingSteps);
};

#endif
