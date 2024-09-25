#include <iostream>
#include <fstream> // Para manejar archivos
#include <cstdlib> // Para rand() y srand()
#include <ctime>   // Para time()
#include <cmath>
#include <vector>  

// Función que genera numeros random.
double random_uniform(double min, double max) {
    return min + static_cast<double>(rand()) / (static_cast<double>(RAND_MAX) / (max - min));
}

// Función Markov-pi
int markov_pi(int N, double delta){
    int n_hits = 0;
    double x = 1.0;
    double y = 1.0;
    for (int i = 0; i < N; ++i){
        double dx = random_uniform(-delta, delta);
        double dy = random_uniform(-delta, delta);
        if (std::abs(x + dx) < 1.0 && std::abs(y + dy) < 1.0){
            x += dx;
            y += dy;
        }
        if ((x * x) + (y * y) < 1.0){
            ++n_hits;
        }
    }
    return n_hits;
}

// Función que calcula la desviación cuadrática media (MSE)
double compute_mse(const std::vector<double>& pi_estimates, double pi_4) {
    double mse = 0.0;
    int num_runs = pi_estimates.size();
    
    for (int i = 0; i < num_runs; ++i) {
        double diff = pi_estimates[i] - pi_4;
        mse += diff * diff;
    }
    
    return mse / num_runs;  // Retorna el promedio de las diferencias al cuadrado
}


int main() {
    srand(time(0)); 
    // Inicializa la semilla para rand()
    std::ofstream output_pi_markov;
    output_pi_markov.open("markov_pi.dat");
    // Con esto puedo probar que se llega a pi/4 con delta = 0.3
    for (int exp = 1; exp <= 8; ++exp){
        int N = static_cast<int>(pow(10,exp));
        for (int i = 0; i < 20; ++i){
            int hits = markov_pi(N, 0.3);
            double pi_estimate = static_cast<double>(hits) / N;
            // std:: cout << pi_estimate << std::endl;
            output_pi_markov << N << " " << pi_estimate << std::endl;
        }
    }
    output_pi_markov.close();

    // std::ofstream output_mse;
    // output_mse.open("mse_markov_pi.dat");

    // const int N_fixed = static_cast<int>(pow(10, 6)); // Usa un N fixed grande
    // const int runs = 20; // Número de simulaciones para cada delta
    // const double pi_4 = M_PI / 4.0;
    // for (double delta = 0.1; delta <= 3.0; delta += 0.01){
    //     std::vector<double> pi_estimates;
    //     for (int i = 0; i < runs; ++i){
    //         int hits = markov_pi(N_fixed, delta);
    //         double pi_estimate = (double)hits / N_fixed;
    //         pi_estimates.push_back(pi_estimate);
    //     }
    //     double mse = compute_mse(pi_estimates, pi_4);
    //     output_mse << delta << " " << mse << std::endl;
    //     std::cout << delta << " " << mse << std::endl;
    // }
    // output_mse.close();

    // std::ofstream rejection_pi_markov;
    // rejection_pi_markov.open("rejection_markov.dat");
    // int N_fixed_2 = static_cast<int>(pow(10,8));
    // for (double delta = 0.0; delta <= 3; delta += 0.01){
    //     int hits = markov_pi(N_fixed_2, delta);
    //     double rejection = 1.0 - (double)hits / N_fixed_2;
    //     double pi_estimate = hits / N_fixed_2;
    //     rejection_pi_markov << delta << " " << rejection << std::endl;
    //     std::cout << delta << " " << rejection << std::endl;
    // }
    // rejection_pi_markov.close();

    
    std::cout << markov_pi(1000, 0.3) << std::endl;
    std::cout << "Todo bien" << std::endl;
    return 0;
}

    // for (double delta = 0.0; delta <= 3.0; delta += 0.01) {
    //     std::vector<double> pi_estimates(runs);
    //     double mse = 0.0;

    //     // Corre num_runs simulaciones para cada delta
    //     for (int i = 0; i < runs; ++i) {
    //         int hits = markov_pi(N_fixed, delta);
    //         pi_estimates[i] = static_cast<double>(hits) / N_fixed;
    //     }

    //     // Calcula la desviación cuadrática media (MSE)
    //     for (int i = 0; i < runs; ++i) {
    //         double diff = pi_estimates[i] - pi_4;
    //         mse += diff * diff;
    //     }
    //     mse /= runs; // Promedio de la desviación cuadrática

    //     // Guarda delta y MSE en el archivo
    //     output_mse << delta << " " << mse << std::endl;
    //     std::cout << "Delta: " << delta << ", MSE: " << mse << std::endl;
    // }

    // output_mse.close();
