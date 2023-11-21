#pragma once
#include "MathFunction.cpp"
#include <random>
#include <vector>

std::random_device random_device;
std::mt19937 random_generator(random_device());

class Gene {
public:
  std::vector<char> chromosomes;
  double fitness = -1;

  Gene(int number_of_chromosomes) {
    std::uniform_int_distribution<int> distribution(0, 1);

    for (int i = 0; i < number_of_chromosomes; i++) {
      int bit = distribution(random_generator);

      if (bit == 0)
        chromosomes.push_back('0');
      else
        chromosomes.push_back('1');
    }
  }

  void flip_chromosome(int index) {
    auto &bit = chromosomes.at(index);

    if (bit == '0')
      bit = '1';
    else
      bit = '0';
  }

  std::vector<double> to_numbers(MathFunction math_function, int dimensions ) {
    std::vector<double> result;
    int chromosomes_per_dimension = chromosomes.size() / dimensions;

    for (int d = 0; d < dimensions; d++) {
      double sum = 0.0;
      for (int i = 0; i < chromosomes_per_dimension; i++) {
        sum = sum * 2 + (chromosomes[d * chromosomes_per_dimension + i] - '0');
      }

      double real_val = math_function.bounds.min +
                        sum * (math_function.bounds.max - math_function.bounds.min) /
                            (pow(2, chromosomes_per_dimension) - 1);
      result.push_back(real_val);
    }

    return result;
  }

  double calculate_function_value(MathFunction math_function, int dimensions) {
    std::vector<double> values = to_numbers(math_function, dimensions);
    return math_function.function(values, dimensions);
  }
};