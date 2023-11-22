#include <algorithm>
#include <iostream>

#include "Chromosome.h"
#include "Population.h"

int main() {
    int dimensions = 1;
    int precision = 5;
    int population_size = 100;
    int generations = 1000;
    double mutation_rate = 0.1;
    double crossover_rate = 0.5;

    auto math_function = rastrigin;

    Population population(population_size, dimensions, precision, math_function);
    Population new_population(population_size, dimensions, precision,
                              math_function);

    int elitism_size =
            static_cast<int>(population_size * 0.05); // 5% pentru elitism

    for (int t = 0; t < generations; t++) {
        population.evaluate_population(math_function);

        for (size_t i = 0; i < population.chromosomes.size(); ++i) {
            for (size_t j = 0; j < population.chromosomes.size() - i - 1; ++j) {
                if (population.chromosomes[j]->fitness < population.chromosomes[j + 1]->fitness) {
                    std::swap(population.chromosomes[j], population.chromosomes[j + 1]);
                }
            }
        }

        // Păstrează cei mai buni 5% indivizi pentru elitism
        new_population.chromosomes.clear();
        new_population.chromosomes.insert(
            new_population.chromosomes.end(), population.chromosomes.begin(),
            population.chromosomes.begin() + elitism_size);

        // Completează restul populației
        auto selected_population =
                new Population(population_size, dimensions, precision, math_function);

        selected_population =
                population.select_population(population_size - elitism_size);
        selected_population->crossover(crossover_rate);
        selected_population->mutation(mutation_rate);

        // Adaugă indivizii selectați, încrucișați și mutați în noua populație
        new_population.chromosomes.insert(new_population.chromosomes.end(),
                                          selected_population->chromosomes.begin(),
                                          selected_population->chromosomes.end());
        // Actualizează populația pentru următoarea iterație
        population = new_population;
    }

    // Găsește cel mai bun individ după ultima iterație
    Chromosome* best_chromosome = nullptr;

    if (!population.chromosomes.empty()) {
        best_chromosome = population.chromosomes[0];
        for (Chromosome* current_chromosome: population.chromosomes) {
            if (current_chromosome->fitness < best_chromosome->fitness) {
                best_chromosome = current_chromosome;
            }
        }
    }

    std::cout << "Best chromosome: ";
    for (char gene: best_chromosome->genes) {
        std::cout << gene;
    }
    std::cout << std::endl;

    std::cout << "Best chromosome value: "
            << best_chromosome->calculate_function_value(math_function,
                                                         dimensions)
            << std::endl;

    std::cout << "Global minimum: " << math_function.global_minimum << std::endl;

    return 0;
}
