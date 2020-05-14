#include "CorrelationCalculator.h"

#include <iomanip>
#include <iostream>

CorrelationCalculator::CorrelationCalculator(std::unique_ptr<TrajectoryFileReader> reader): nTimeframes{reader->getNTimeframes()}, velocities{reader->getVelocities()}, fReader{std::move(reader)}
{
    nAtoms = velocities.at(0).size();
    correlations.resize(nTimeframes);
}

void CorrelationCalculator::printCorrelation()
{
    calculateCorrelation();
    for (std::size_t i = 0; i < correlations.size(); ++i) {
        std::cout << i << std::setw(10) << correlations.at(i) << std::endl;
    }
}

void CorrelationCalculator::calculateCorrelation()
{
    if (fReader->getCalculationType() == 0) {
        for (std::size_t i = 0; i < nTimeframes; ++i) {
            for (std::size_t j = 0; j < nAtoms; ++j) {
                double sum{0};
                for (std::size_t k = 0; k < nTimeframes-i; ++k) {
                    sum += velocities.at(i).at(k)*velocities.at(i).at(k+i);
                }
                correlations[i] += sum/(nTimeframes-i);
            }
            correlations[i] /= nAtoms;
        }
    } else {
        std::cout << "Unimplemented for the time being.\n";
    }
}
