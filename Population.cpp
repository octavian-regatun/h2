#pragma once
#include "Chromosome.cpp"
#include "MathFunction.cpp"
#include <vector>

class Population{
public:
  std::vector<Chromosome *> chromosomes;
  int dimensions;
  int precision;
  MathFunction math_function;

  Population(int number_of_chromosomes, int dimensions, int precision, MathFunction math_function)
      : dimensions(dimensions), precision(precision), math_function(math_function) {
    int number_of_genes = calculate_length();

    for (int i = 0; i < number_of_chromosomes; i++) {
      genes.push_back(new Chromosome(number_of_genes));
    }
  }

  Population(const Population& other)
    : dimensions(other.dimensions), precision(other.precision), math_function(other.math_function) {
    chromosomes.reserve(other.chromosomes.size());
    for (auto& chromosome : other.chromosomes) {
        chromosomes.push_back(new Chromosome(*chromosome));
    }
  }


  int calculate_length() {
    return static_cast<int>(std::ceil(
        dimensions * std::log2(pow(10, precision) *
                               (math_function.bounds.max - math_function.bounds.min))));
  }

  double evaluate_population(MathFunction math_function) {
    double C = 1000.0; 

    for (auto& individual : chromosomes) {
        double function_value = individual.calculate_function_value(math_function, dimensions);
        if (function_value < 0) {
            individual.fitness = function_value + C; 
        } else {
            individual.fitness = 1.0 / (function_value + 1.0); 
        }
    }
  }

  void mutation(double mutation_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    for (auto& individual : chromosomes) {
        for (int index = 0; index < dimensions; index++) {
            if (dis(gen) < mutation_rate) {
                individual.flip_gene(index); // Flip bit
            }
        }
    }
  }

  void crossover_population(double crossover_rate) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<> point_dis(1, chromosomes[0].genes.size() - 2); // Exclude capetele

    for (size_t i = 0; i < chromosomes.size() - 1; i += 2) {
        if (dis(gen) < crossover_rate) {
            int crossover_point = point_dis(gen);
            for (int j = crossover_point; j < chromosomes[i].genes.size(); j++) {
                std::swap(chromosomes[i].genes[j], chromosomes[i + 1].genes[j]);
            }
        }
    }
  }

  Population& select_population(int num_selected) {
    double total_fitness = std::accumulate(chromosomes.begin(), chromosomes.end(), 0.0, [](double sum, const Individual& ind) {
        return sum + ind.fitness;
    });

    std::vector<double> probabilities(chromosomes.size());
    for (size_t i = 0; i < chromosomes.size(); i++) {
        probabilities[i] = chromosomes[i].fitness / total_fitness;
    }

    std::vector<double> cumulative_probabilities(chromosomes.size() + 1, 0.0);
    std::partial_sum(probabilities.begin(), probabilities.end(), cumulative_probabilities.begin() + 1);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<Chromosome *> selected_population;
    for (int i = 0; i < num_selected; i++) {
        double r = dis(gen);
        auto it = std::upper_bound(cumulative_probabilities.begin(), cumulative_probabilities.end(), r);
        int index = std::distance(cumulative_probabilities.begin(), it) - 1;

        selected_population.push_back(chromosomes[index]);
    }

    Population new_population(num_selected, dimensions, precision, math_function);
    for (int i = 0; i < num_selected; ++i) {
        new_population.chromosomes[i] = selected_population[i];
    }

    return new_population;
  }

};