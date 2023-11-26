#include <algorithm>
#include <iostream>

#include "Chromosome.h"
#include "Population.h"

int main() {
    // ----------CONSTANTS----------
    int dimensions = 30;
    int precision = 5;
    int population_size = 200;
    int generations = 2000;
    double mutation_rate = 1.0 / 800.0;
    double crossover_rate = 0.6;
    double elitism_rate = 0.05;
    MathFunction math_function = rastrigin;

    // ----------CREATE POPULATION----------
    Population population(population_size, dimensions, precision, rastrigin);

    for (int i = 1; i <= generations; i++) {
        population.normalize_fitness();
        population.calculate_cumsum();
        population.crossover(crossover_rate, elitism_rate, mutation_rate);
        std::cout << "Generation " << i << " done!" << std::endl;
    }


    // ----------DEBUG----------
    population.normalize_fitness();
    population.calculate_cumsum();
    int index = 0;
    for (auto&chromosome: population.chromosomes) {
        printf("Chromosome %d:\n", ++index);
        printf("Function Value: %.5f\n", chromosome.calculate_function_value(math_function, dimensions));
        printf("Fitness: %.3f\n", chromosome.fitness);
        printf("Cumsum: %.3f\n", chromosome.cumsum);
        printf("\n\n");
    }

    // display the best chromosome
    auto best_chromosome = population.chromosomes[0];
    printf("Function Global Minimum: %.5f\n", math_function.global_minimum);
    printf("Best Chromosome:\n");
    printf("Function Value: %.5f\n", best_chromosome.calculate_function_value(math_function, dimensions));
    printf("Fitness: %.3f\n", best_chromosome.fitness);
    printf("Cumsum: %.3f\n", best_chromosome.cumsum);
    printf("\n\n");

    printf("Program Done Running!");

    return 0;
}
