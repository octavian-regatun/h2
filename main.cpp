#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <vector>

double PI = 3.14159265358979323846;

struct bounds {
  double min, max;
};

double rastrigin(const std::vector<double> &x, int dimensions) {
  double sum = 10 * dimensions;
  for (int i = 0; i < dimensions; i++)
    sum += x[i] * x[i] - 10 * cos(2 * PI * x[i]);
  return sum;
}

double michalewicz(const std::vector<double> &x, int dimensions) {
  double sum = 0;
  for (int i = 0; i < dimensions; i++)
    sum -= sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / PI), 2 * 10);
  return sum;
}

double dejong1(const std::vector<double> &x, int dimensions) {
  double sum = 0.0;
  for (int i = 0; i < dimensions; i++)
    sum += x[i] * x[i];
  return sum;
}

double schwefel(const std::vector<double> &x, int dimensions) {
  double sum = 0.0;
  for (int i = 0; i < dimensions; i++) {
    sum += x[i] * sin(sqrt(std::abs(x[i])));
  }
  return 418.9829 * dimensions - sum;
}

int calculate_length(const bounds &bounds, int dimensions, int precision) {
  return static_cast<int>(std::ceil(
      dimensions * std::log2(pow(10, precision) * (bounds.max - bounds.min))));
}

std::random_device random_device;
std::mt19937 random_generator(random_device());

std::vector<char> generate_random_bits(int length) {
  std::uniform_int_distribution<int> distribution(0, 1);

  std::vector<char> bits;

  for (int i = 0; i < length; i++) {
    int bit = distribution(random_generator);

    if (bit == 0)
      bits.push_back('0');
    else
      bits.push_back('1');
  }

  return bits;
}

void flipBit(std::vector<char> &bits, int index) {
  auto &bit = bits.at(index);

  if (bit == '0')
    bit = '1';
  else
    bit = '0';
}

std::vector<double> bits_to_number(const std::vector<char> &bits,
                                   const bounds &bounds, int length,
                                   int dimensions) {
  std::vector<double> result;
  int bits_per_dimension = length / dimensions;

  for (int d = 0; d < dimensions; d++) {
    double sum = 0.0;
    for (int i = 0; i < bits_per_dimension; i++) {
      sum = sum * 2 + (bits[d * bits_per_dimension + i] - '0');
    }

    double real_val = bounds.min + sum * (bounds.max - bounds.min) /
                                       (pow(2, bits_per_dimension) - 1);
    result.push_back(real_val);
  }

  return result;
}

double calculate_function(const std::vector<char> &bits,
                          double (*func)(const std::vector<double> &, int),
                          const bounds &bounds, int dimensions, int length) {
  std::vector<double> values = bits_to_number(bits, bounds, length, dimensions);
  return func(values, dimensions);
}

// Structura pentru a reține un individ și scorul său
struct Individual {
    std::vector<char> bits;
    double fitness;
};

// Funcția de generare a populației inițiale
std::vector<Individual> generate_initial_population(int population_size, int length) {
    std::vector<Individual> population;
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < population_size; i++) {
        std::vector<char> chromosome(length);
        std::uniform_int_distribution<> dis(0, 1);
        std::generate(chromosome.begin(), chromosome.end(), [&]() { return dis(gen); });

        population.push_back({chromosome, 0.0}); // Fitness inițial 0.0
    }

    return population;
}


void evaluate_population(std::vector<Individual>& population, double (*func)(const std::vector<double>&, int), const bounds& bounds, int dimensions) {
    double C = 1000.0; 

    for (auto& individual : population) {
        double function_value = calculate_function(individual.bits, func, bounds, dimensions);
        if (function_value < 0) {
            individual.fitness = function_value + C; 
        } else {
            individual.fitness = 1.0 / (function_value + 1.0); 
        }
    }
}

//Roata norocului
std::vector<Individual> select_population(const std::vector<Individual>& population, int num_selected) {
    double total_fitness = std::accumulate(population.begin(), population.end(), 0.0, [](double sum, const Individual& ind) {
        return sum + ind.fitness;
    });

    std::vector<double> probabilities(population.size());
    for (size_t i = 0; i < population.size(); i++) {
        probabilities[i] = population[i].fitness / total_fitness;
    }

    std::vector<double> cumulative_probabilities(population.size() + 1, 0.0);
    std::partial_sum(probabilities.begin(), probabilities.end(), cumulative_probabilities.begin() + 1);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<Individual> selected_population;
    for (int i = 0; i < num_selected; i++) {
        double r = dis(gen);
        auto it = std::upper_bound(cumulative_probabilities.begin(), cumulative_probabilities.end(), r);
        int index = std::distance(cumulative_probabilities.begin(), it) - 1;

        selected_population.push_back(population[index]);
    }

    return selected_population;
}

// mutation
void mutate_population(std::vector<Individual>& population, double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (auto& individual : population) {
        for (auto& bit : individual.bits) {
            if (dis(gen) < mutation_rate) {
                bit = bit ^ 1; // Flip bit
            }
        }
    }
}

// crossover
void crossover_population(std::vector<Individual>& population, double crossover_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> point_dis(1, population[0].bits.size() - 2); // Exclude capetele

    for (size_t i = 0; i < population.size() - 1; i += 2) {
        if (dis(gen) < crossover_rate) {
            int crossover_point = point_dis(gen);
            for (int j = crossover_point; j < population[i].bits.size(); j++) {
                std::swap(population[i].bits[j], population[i + 1].bits[j]);
            }
        }
    }
}

// Algoritm genetic crossover plus elitism
std::vector<char> genetic_algorithm(double (*func)(const std::vector<double>&, int), const bounds& bounds, int dimensions, int population_size, int iterations, double mutation_rate, double crossover_rate) {
    std::vector<Individual> population = generate_initial_population(population_size, dimensions);
    std::vector<Individual> new_population;
    int elitism_size = static_cast<int>(population_size * 0.05); // 5% pentru elitism

    for (int t = 0; t < iterations; t++) {
        evaluate_population(population, func, bounds, dimensions);

        
        std::sort(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
            return a.fitness > b.fitness; // Sortează descrescător pentru minimizare
        });

        // Păstrează cei mai buni 5% indivizi pentru elitism
        new_population.clear();
        new_population.insert(new_population.end(), population.begin(), population.begin() + elitism_size);

        // Completează restul populației
        std::vector<Individual> selected_population = select_population(population, population_size - elitism_size);
        crossover_population(selected_population, crossover_rate);
        mutate_population(selected_population, mutation_rate);

        // Adaugă indivizii selectați, încrucișați și mutați în noua populație
        new_population.insert(new_population.end(), selected_population.begin(), selected_population.end());

        // Actualizează populația pentru următoarea iterație
        population = new_population;
    }

    // Găsește cel mai bun individ după ultima iterație
    auto best_it = std::min_element(population.begin(), population.end(), [](const Individual& a, const Individual& b) {
        return a.fitness < b.fitness;
    });

    return best_it->bits;
}




int main() {
  bounds dejong1_bounds = {-5.12, 5.12};
  bounds schwefel_bounds = {-500, 500};
  bounds rastrigin_bounds = {-5.12, 5.12};
  bounds michalewicz_bounds = {0, M_PI};

  
}