#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <omp.h>
#include <random>
#include <string>
#include <vector>

double PI = 3.14159265358979323846;
int LENGTH, PRECISION = 5;

struct bounds {
  double min, max;
};

double rastrigin(const std::vector<double> &x, int dimensions) {
  double sum = 10 * dimensions;
  for (int i = 0; i < dimensions; i++)
    sum += x[i] * x[i] - 10 * cos(2 * PI * x[i]);
  return sum;
}

double michalewicz(const std::vector<double> &x, int dimensions) {
  double sum = 0;
  for (int i = 0; i < dimensions; i++)
    sum -= sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / PI), 2 * 10);
  return sum;
}

double dejong1(const std::vector<double> &x, int dimensions) {
  double sum = 0.0;
  for (int i = 0; i < dimensions; i++)
    sum += x[i] * x[i];
  return sum;
}

double schwefel(const std::vector<double> &x, int dimensions) {
  double sum = 0.0;
  for (int i = 0; i < dimensions; i++) {
    sum += x[i] * sin(sqrt(abs(x[i])));
  }
  return 418.9829 * dimensions - sum;
}

int calculate_length(const bounds &bounds, int dimensions, int precision) {
  return static_cast<int>(std::ceil(
      dimensions * std::log2(pow(10, precision) * (bounds.max - bounds.min))));
}

std::random_device random_device;
std::mt19937 random_generator(omp_get_thread_num() + random_device());

std::vector<char> generate_random_bits(int length) {
  std::uniform_int_distribution<int> distribution(0, 1);

  std::vector<char> bits;

  for (int i = 0; i < length; i++) {
    int bit = distribution(random_generator);

    if (bit == 0)
      bits.push_back('0');
    else
      bits.push_back('1');
  }

  return bits;
}

void flipBit(std::vector<char> &bits, int index) {
  auto &bit = bits.at(index);

  if (bit == '0')
    bit = '1';
  else
    bit = '0';
}

std::vector<double> bits_to_number(const std::vector<char> &bits,
                                   const bounds &bounds, int length,
                                   int dimensions) {
  std::vector<double> result;
  int bits_per_dimension = length / dimensions;

  for (int d = 0; d < dimensions; d++) {
    double sum = 0.0;
    for (int i = 0; i < bits_per_dimension; i++) {
      sum = sum * 2 + (bits[d * bits_per_dimension + i] - '0');
    }

    double real_val = bounds.min + sum * (bounds.max - bounds.min) /
                                       (pow(2, bits_per_dimension) - 1);
    result.push_back(real_val);
  }

  return result;
}

double calculate_function(const std::vector<char> &bits,
                          double (*func)(const std::vector<double> &, int),
                          const bounds &bounds, int dimensions, int length) {
  std::vector<double> values = bits_to_number(bits, bounds, length, dimensions);
  return func(values, dimensions);
}

int main() { std::cout << "Hello world"; }