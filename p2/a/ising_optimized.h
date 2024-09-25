#include <bitset>
#include <iostream>

class IsingModelOptimized {
private:
    int N;  // Tamaño de la red
    long long spins;  // Representación de la red de spins con bits

public:
    IsingModelOptimized(int N) : N(N) {}

    // Convierte un entero a configuración de spins
    std::bitset<36> toBits(long long config) {
        return std::bitset<36>(config);
    }

    // Calcula la energía de una configuración usando bits
    int energyFromBits(long long config, int J) {
        std::bitset<36> spins = toBits(config);
        int energy = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                int current = i * N + j;
                int right = i * N + ((j + 1) % N);  // Vecino a la derecha con PBC
                int down = ((i + 1) % N) * N + j;  // Vecino de abajo con PBC

                energy -= J * (spins[current] == spins[right] ? 1 : -1);
                energy -= J * (spins[current] == spins[down] ? 1 : -1);
            }
        }
        return energy;
    }
};
