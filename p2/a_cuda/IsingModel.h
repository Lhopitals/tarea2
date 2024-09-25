#ifndef ISINGMODEL_H
#define ISINGMODEL_H

#include <vector>

class IsingModel {
private:
    int N; // Tamaño de la red
    std::vector<std::vector<int>> spins; // Guardar spins en CPU

public:
    // Constructor de la clase
    IsingModel(int N);

    // Destructor de la clase
    ~IsingModel();

    // Método para enumerar todas las configuraciones posibles de spins
    std::vector<std::vector<std::vector<int>>> EnumerateIsing() const;

    // Función para calcular la energía con condiciones periódicas en la GPU
    void calculateEnergiesWithEnumerateIsingCUDA(int J, std::vector<int>& energiesNonPeriodic, std::vector<int>& energiesPeriodic);
};

#endif
