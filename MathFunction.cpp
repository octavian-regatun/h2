#pragma once
#include <math.h>
#include <vector>

double PI = 3.14159265358979323846;

struct Bounds {
  double min, max;
};

class MathFunction {
public:
  double (*function)(const std::vector<double> &, int);
  Bounds bounds;
  double global_minimum;

  MathFunction(double (*f)(const std::vector<double> &, int), Bounds bounds,
           double global_minimum)
      : function(f), bounds(bounds), global_minimum(global_minimum) {}
};

auto rastrigin = MathFunction(
    [](const std::vector<double> &x, int dimensions) {
      double sum = 10 * dimensions;
      for (int i = 0; i < dimensions; i++)
        sum += x[i] * x[i] - 10 * cos(2 * PI * x[i]);
      return sum;
    },
    {-5.12, 5.12}, 0);

auto michalewicz = MathFunction(
    [](const std::vector<double> &x, int dimensions) {
      double sum = 0;
      for (int i = 0; i < dimensions; i++)
        sum -= sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / PI), 2 * 10);
      return sum;
    },
    {0, PI}, -4.687658);

auto dejong1 = MathFunction(
    [](const std::vector<double> &x, int dimensions) {
      double sum = 0.0;
      for (int i = 0; i < dimensions; i++)
        sum += x[i] * x[i];
      return sum;
    },
    {-5.12, 5.12}, 0);

auto schwefel = MathFunction(
    [](const std::vector<double> &x, int dimensions) {
      double sum = 0.0;
      for (int i = 0; i < dimensions; i++) {
        sum += x[i] * sin(sqrt(std::abs(x[i])));
      }
      return 418.9829 * dimensions - sum;
    },
    {-500, 500}, 0);
