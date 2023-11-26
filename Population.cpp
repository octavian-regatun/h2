//
// Created by Octavian Regatun on 22-Nov-23.
//

#include "Population.h"
#include "Chromosome.h"
#include "MathFunction.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <vector>
#include <random>

Population::Population(int number_of_chromosomes, int dimensions, int precision,
                       MathFunction math_function): math_function(math_function), dimensions(dimensions),
                                                    precision(precision) {
    int number_of_genes = calculate_length();

    for (int i = 0; i < number_of_chromosomes; i++) {
        auto chromosome = new Chromosome(number_of_genes);
        chromosomes.push_back(*chromosome);
    }
}

int Population::calculate_length() {
    return static_cast<int>(
        std::ceil(dimensions *
                  std::log2(std::pow(10, precision) * (math_function.bounds.max -
                                                       math_function.bounds.min))));
}

void Population::normalize_fitness() {
    for (auto&chromosome: chromosomes) {
        chromosome.calculate_fitness(math_function, dimensions);
        chromosome.fitness = 1 / chromosome.fitness;
    }

    double fitness_sum = 0;
    for (auto&chromosome: chromosomes) {
        fitness_sum += chromosome.fitness;
    }

    for (auto&chromosome: chromosomes) {
        chromosome.fitness /= fitness_sum;
    }
}

void Population::calculate_cumsum() {
    // sort them based on fitness value (descending order)
    std::sort(chromosomes.begin(), chromosomes.end(),
              [](const Chromosome&lhs, const Chromosome&rhs) {
                  return lhs.fitness > rhs.fitness;
              });

    for (int i = 0; i < chromosomes.size(); i++) {
        chromosomes[i].cumsum = 0;
        for (int j = i; j < chromosomes.size(); j++) {
            chromosomes[i].cumsum += chromosomes[j].fitness;
        }
    }
}

void Population::crossover(double crossover_rate, double elitism_rate, double mutation_rate) {
    auto elites_chromosomes = elitism(elitism_rate);
    std::vector<Chromosome> new_chromosomes = elites_chromosomes;

    for (int i = 0; i < (chromosomes.size() - elites_chromosomes.size()) / 2; i++) {
        std::uniform_real_distribution<double> distribution(1, calculate_length() - 1);
        auto random_generator = std::mt19937(std::random_device()());
        int crossover_point = static_cast<int>(distribution(random_generator));

        auto parent1 = pick_parent();
        // TODO: write a more efficient way to pick parent2
        Chromosome* parent2 = nullptr;
        do {
            parent2 = pick_parent();
        }
        while (parent2 == nullptr && parent2->fitness == parent1->fitness);

        auto childs = Chromosome::crossover(parent1, parent2, crossover_point, crossover_rate);

        childs[0].mutate(mutation_rate);
        childs[1].mutate(mutation_rate);

        new_chromosomes.insert(new_chromosomes.end(), childs.begin(), childs.end());
    }

    chromosomes = new_chromosomes;
}

Chromosome* Population::pick_parent() {
    std::uniform_real_distribution<double> distribution(0, 1);
    auto random_generator = std::mt19937(std::random_device()());
    double random = distribution(random_generator);

    for (auto&chromosome: chromosomes) {
        if (random >= chromosome.cumsum) {
            return &chromosome;
        }
    }
    return &chromosomes[chromosomes.size() - 1];
}

void Population::mutate(double mutation_rate) {
    for (auto&chromosome: chromosomes) {
        chromosome.mutate(mutation_rate);
    }
}

std::vector<Chromosome> Population::elitism(double elitism_rate) {
    int number_of_elites = floor(chromosomes.size() * elitism_rate);
    std::vector<Chromosome> elite_chromosomes;

    for (int i = 0; i < number_of_elites; i++) {
        elite_chromosomes.push_back(chromosomes[i]);
    }

    return elite_chromosomes;
}
