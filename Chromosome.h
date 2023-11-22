#ifndef CHROMOSOME_H
#define CHROMOSOME_H
#include <vector>
#include "MathFunction.h"


class Chromosome {
public:
    std::vector<char> genes;
    double fitness = -1;

    explicit Chromosome(int number_of_genes);

    void flip_gene(int index);

    std::vector<double> to_numbers(MathFunction math_function, int dimensions);

    double calculate_function_value(MathFunction math_function, int dimensions);
};


#endif //CHROMOSOME_H