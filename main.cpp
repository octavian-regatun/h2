#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "Chromosome.cpp"
#include "MathFunction.cpp"
#include "Population.cpp"

double PI = 3.14159265358979323846;

std::random_device random_device;
std::mt19937 random_generator(random_device());

#define DIMENSIONS 30
#define PRECISION 5
#define POPULATION_SIZE 100
#define GENERATIONS 1000
#define MUTATION_RATE 0.01
#define CROSSOVER_RATE 0.5
#define NUMBER_OF_GENES 10

// Algoritm genetic crossover plus elitism
std::vector<char> genetic_algorithm(MathFunction math_function, int dimensions, int population_size, int iterations, double mutation_rate, double crossover_rate) {
    Population population(population_size, dimensions, PRECISION, math_function);
    Population new_population(population_size, dimensions, PRECISION, math_function);
    int elitism_size = static_cast<int>(population_size * 0.05); // 5% pentru elitism

    for (int t = 0; t < iterations; t++) {
        population.evaluate_population(math_function);

        std::sort(population.chromosomes.begin(), population.chromosomes.end(), [](const Individual& a, const Individual& b) {
            return a.fitness > b.fitness; // Sortează descrescător pentru minimizare
        });

        // Păstrează cei mai buni 5% indivizi pentru elitism
        new_population.chromosomes.clear();
        new_population.chromosomes.insert(new_population.chromosomes.end(), population.chromosomes.begin(), population.chromosomes.begin() + elitism_size);

        // Completează restul populației
        Population aux(population_size, dimensions, PRECISION, math_function);
        aux.select_population(population_size - elitism_size);
        aux.crossover(crossover_rate);
        aux.mutation(mutation_rate);

        // Adaugă indivizii selectați, încrucișați și mutați în noua populație
        new_population.chromosomes.insert(new_population.chromosomes.end(), selected_population.chromosomes.begin(), selected_population.chromosomes.end());
        // Actualizează populația pentru următoarea iterație
        population = new_population;
    }

    // Găsește cel mai bun individ după ultima iterație
    auto best_it = std::min_element(population.chromosomes.begin(), population.chromosomes.end(), [](const Individual& a, const Individual& b) {
        return a.fitness < b.fitness;
    });

    return best_it->genes;
}



