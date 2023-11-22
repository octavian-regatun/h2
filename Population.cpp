//
// Created by Octavian Regatun on 22-Nov-23.
//

#include "Population.h"
#include "Chromosome.h"
#include "MathFunction.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <random>

Population::Population(int number_of_chromosomes, int dimensions, int precision,
                       MathFunction math_function): math_function(math_function), dimensions(dimensions),
                                                    precision(precision) {
    int number_of_genes = calculate_length();

    for (int i = 0; i < number_of_chromosomes; i++) {
        chromosomes.push_back(new Chromosome(number_of_genes));
    }
}

int Population::calculate_length() {
    return static_cast<int>(
        std::ceil(dimensions *
                  std::log2(std::pow(10, precision) * (math_function.bounds.max -
                                                       math_function.bounds.min))));
}

void Population::evaluate_population(MathFunction math_function) {
    double C = 1000.0;

    for (auto&chromosome: chromosomes) {
        double function_value =
                chromosome->calculate_function_value(math_function, dimensions);
        if (function_value < 0) {
            chromosome->fitness = function_value + C;
        }
        else {
            chromosome->fitness = 1.0 / (function_value + 1.0);
        }
    }
}

void Population::mutation(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (auto&chromosome: chromosomes) {
        for (int index = 0; index < dimensions; index++) {
            if (dis(gen) < mutation_rate) {
                chromosome->flip_gene(index);
            }
        }
    }
}

void Population::crossover(double crossover_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> point_dis(1, chromosomes[0]->genes.size() -
                                                 2); // Exclude capetele

    for (size_t i = 0; i < chromosomes.size() - 1; i += 2) {
        if (dis(gen) < crossover_rate) {
            int crossover_point = point_dis(gen);
            for (int j = crossover_point; j < chromosomes[i]->genes.size(); j++) {
                std::swap(chromosomes[i]->genes[j], chromosomes[i + 1]->genes[j]);
            }
        }
    }
}

Population* Population::select_population(int num_selected) {
    // Calculul fitness-ului total
    double total_fitness =
            std::accumulate(chromosomes.begin(), chromosomes.end(), 0.0,
                            [](double sum, Chromosome* chromosome) {
                                return sum + chromosome->fitness;
                            });

    // Calculul probabilităților pentru fiecare cromozom
    std::vector<double> probabilities(chromosomes.size());
    for (size_t i = 0; i < chromosomes.size(); i++) {
        probabilities[i] = chromosomes[i]->fitness / total_fitness;
    }

    // Calculul probabilităților cumulative
    std::vector<double> cumulative_probabilities(chromosomes.size() + 1, 0.0);
    std::partial_sum(probabilities.begin(), probabilities.end(),
                     cumulative_probabilities.begin() + 1);

    // Inițializarea generatorului de numere aleatorii
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    // Selecția cromozomilor pe baza probabilităților cumulative
    std::vector<Chromosome *> selected_population;
    for (int i = 0; i < num_selected; i++) {
        double r = dis(gen);
        auto it = std::upper_bound(cumulative_probabilities.begin(),
                                   cumulative_probabilities.end(), r);
        int index = std::distance(cumulative_probabilities.begin(), it) - 1;

        selected_population.push_back(chromosomes[index]);
    }

    // Crearea unei noi populații cu cromozomii selectați
    Population* new_population = new Population(num_selected, dimensions, precision, math_function);
    for (int i = 0; i < num_selected; ++i) {
        new_population->chromosomes[i] = new Chromosome(*selected_population[i]);
    }

    return new_population;
}
