#ifndef CORRELATIONCALCULATOR_H
#define CORRELATIONCALCULATOR_H

#include "TrajectoryFileReader.h"

#include <memory>
#include <vector>

class CorrelationCalculator
{
public:
    CorrelationCalculator(std::unique_ptr<TrajectoryFileReader> reader);

    void printCorrelation();
private:
    std::size_t nAtoms;
    std::size_t nTimeframes;
    std::vector<double> correlations;
    std::vector<std::vector<double>> velocities;
    std::unique_ptr<TrajectoryFileReader> fReader;

    void calculateCorrelation();
};

#endif // CORRELATIONCALCULATOR_H
