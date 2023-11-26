//
// Created by Octavian Regatun on 22-Nov-23.
//

#ifndef POPULATION_H
#define POPULATION_H
#include <vector>
#include "Chromosome.h"

class Population {
public:
    std::vector<Chromosome> chromosomes;
    int dimensions;
    int precision;
    MathFunction math_function;

    Population(int number_of_chromosomes, int dimensions, int precision, MathFunction math_function);

    int calculate_length();

    void normalize_fitness();

    void calculate_cumsum();

    Chromosome* pick_parent();

    void crossover(double crossover_rate, double elitism_rate,double mutation_rate);

    void mutate(double mutation_rate);

    std::vector<Chromosome> elitism(double elitism_rate);
};

#endif //POPULATION_H
