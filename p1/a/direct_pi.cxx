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

// Función que calcula el valor de pi.
int direct_pi(int N) {
    int n_hits = 0;
    for (int i = 0; i < N; ++i) {
        double x = random_uniform(-1.0, 1.0);
        double y = random_uniform(-1.0, 1.0);
        if (( (x * x) + (y * y ) ) < 1.0){
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
    srand(time(0)); // Inicializa la semilla para rand()
    std::ofstream output_pi;
    output_pi.open("direct_pi.dat");
    for (int exp = 1; exp <= 8; ++exp) {
        int N = static_cast<int>(pow(10, exp));
        for (int i = 0; i < 20; ++i) {
            int hits = direct_pi(N);
            double pi_estimate =  (double)hits / N; // pi/4
            // double pi_estimate = 4.0 * hits / N; // Acá calculo pi
            output_pi << N << " " << pi_estimate << std::endl;
        }
    }
    output_pi.close();

    std::ofstream output_mse;
    output_mse.open("mse_direct_pi_1.dat");
    for (int exp = 1; exp <= 8; ++exp) {
        int N_fixed = static_cast<int>(pow(10, exp));
        const int runs = 20;
        const double pi_4 = M_PI / 4.0;
        std::vector<double> pi_estimates;
        for (int i = 0; i < runs; ++i) {
            int hits = direct_pi(N_fixed);
            double pi_estimate = (double)hits / N_fixed;
            pi_estimates.push_back(pi_estimate);
        double mse = compute_mse(pi_estimates, pi_4);
        output_mse << N_fixed << " " << mse << std::endl;
        std::cout << N_fixed << " " << mse << std::endl;
        }
    }
    return 0;
}