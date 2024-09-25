#ifndef ISINGMODEL_H
#define ISINGMODEL_H

#include <vector>

class IsingModel {
private:
    int N; // Tamaño de la red
    std::vector<std::vector<int>> spins;

public:
    // Constructor de la clase
    IsingModel(int N);

    // Función para crear una red de spins aleatorios (+1 o -1)
    std::vector<std::vector<int>> createSpinsGrid();

    // Función para imprimir la configuración de spins en formato de matriz
    void printSpins() const;

    // Función para calcular la energía con condiciones de borde periódicas
    int energyIsingPeriodic(int J, const std::vector<std::vector<int> >& spins) const;

    // Función para calcular la energía con condiciones de borde no periódicas
    int energyIsingNonPeriodic(int J, const std::vector<std::vector<int> >& spins) const;

    // Función Enumerate-Ising para generar todas las configuraciones
    std::vector<std::vector<std::vector<int>>> EnumerateIsing() const;

    // Función gray-flip para cambiar un solo spin
    void gray_flip(const std::vector<int>& previousGray, const std::vector<int>& currentGray);
};

#endif
