//
// Created by Octavian Regatun on 22-Nov-23.
//

#ifndef MATHFUNCTION_H
#define MATHFUNCTION_H
#include <vector>

struct Bounds {
    double min, max;
};

class MathFunction {
public:
    double (*function)(const std::vector<double>&, int) = nullptr;

    Bounds bounds = {};
    double global_minimum = 0;

    MathFunction(double (*function)(const std::vector<double>&, int), Bounds bounds, double global_minimum);
};

extern MathFunction rastrigin;
extern MathFunction michalewicz;
extern MathFunction dejong1;
extern MathFunction schwefel;

#endif //MATHFUNCTION_H
