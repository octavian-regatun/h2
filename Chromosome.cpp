//
// Created by Octavian Regatun on 22-Nov-23.
//

#include "Chromosome.h"
#include "MathFunction.h"
#include <random>
#include <vector>

std::random_device random_device;
std::mt19937 random_generator(random_device());

Chromosome::Chromosome(int number_of_genes) {
    std::uniform_int_distribution<int> distribution(0, 1);

    for (int i = 0; i < number_of_genes; i++) {
        int bit = distribution(random_generator);

        if (bit == 0)
            genes.push_back('0');
        else
            genes.push_back('1');
    }
};

void Chromosome::flip_gene(int index) {
    auto&bit = genes.at(index);

    if (bit == '0')
        bit = '1';
    else
        bit = '0';
};

std::vector<double> Chromosome::to_numbers(MathFunction math_function, int dimensions) {
    std::vector<double> result;
    int genes_per_dimension = genes.size() / dimensions;

    for (int d = 0; d < dimensions; d++) {
        double sum = 0.0;
        for (int i = 0; i < genes_per_dimension; i++) {
            sum = sum * 2 + (genes[d * genes_per_dimension + i] - '0');
        }

        double real_val = math_function.bounds.min +
                          sum * (math_function.bounds.max - math_function.bounds.min) /
                          (pow(2, genes_per_dimension) - 1);
        result.push_back(real_val);
    }

    return result;
};

double Chromosome::calculate_function_value(MathFunction math_function, int dimensions) {
    std::vector<double> values = to_numbers(math_function, dimensions);
    return math_function.function(values, dimensions);
};
