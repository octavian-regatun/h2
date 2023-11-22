//
// Created by Octavian Regatun on 22-Nov-23.
//

#ifndef POPULATION_H
#define POPULATION_H
#include <vector>
#include "Chromosome.h"

class Population {
public:
    std::vector<Chromosome *> chromosomes;
    int dimensions{};
    int precision{};
    MathFunction math_function;

    Population(int number_of_chromosomes, int dimensions, int precision, MathFunction math_function);

    int calculate_length();

    void evaluate_population(MathFunction math_function);

    void mutation(double mutation_rate);

    void crossover(double crossover_rate);

    Population* select_population(int num_selected);
};

#endif //POPULATION_H
