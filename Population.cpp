#pragma once
#include "Gene.cpp"
#include <vector>

class Population {
public:
  std::vector<Gene *> genes;
  int dimensions;
  int precision;
  MathFunction math_function;

  Population(int size, int dimensions, int precision, MathFunction math_function)
      : dimensions(dimensions), precision(precision), math_function(math_function) {
    int number_of_chromosomes = calculate_length();

    for (int i = 0; i < size; i++) {
      genes.push_back(new Gene(number_of_chromosomes));
    }
  }

  int calculate_length() {
    return static_cast<int>(std::ceil(
        dimensions * std::log2(pow(10, precision) *
                               (math_function.bounds.max - math_function.bounds.min))));
  }
};