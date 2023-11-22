#include "MathFunction.h"
#include <cmath>
#include <vector>

double PI = 3.14159265358979323846;

MathFunction::MathFunction(double (*function)(const std::vector<double>&, int), Bounds bounds,
                           double global_minimum): function(function), bounds(bounds), global_minimum(global_minimum) {
};

MathFunction rastrigin = MathFunction(
    [](const std::vector<double>&x, int dimensions) {
        double sum = 10 * dimensions;
        for (int i = 0; i < dimensions; i++)
            sum += x[i] * x[i] - 10 * cos(2 * PI * x[i]);
        return sum;
    },
    {-5.12, 5.12}, 0);

MathFunction michalewicz = MathFunction(
    [](const std::vector<double>&x, int dimensions) {
        double sum = 0;
        for (int i = 0; i < dimensions; i++)
            sum -= sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / PI), 2 * 10);
        return sum;
    },
    {0, PI}, -4.687658);

MathFunction dejong1 = MathFunction(
    [](const std::vector<double>&x, int dimensions) {
        double sum = 0.0;
        for (int i = 0; i < dimensions; i++)
            sum += x[i] * x[i];
        return sum;
    },
    {-5.12, 5.12}, 0);

MathFunction schwefel = MathFunction(
    [](const std::vector<double>&x, int dimensions) {
        double sum = 0.0;
        for (int i = 0; i < dimensions; i++) {
            sum += x[i] * sin(sqrt(std::abs(x[i])));
        }
        return 418.9829 * dimensions - sum;
    },
    {-500, 500}, 0);
