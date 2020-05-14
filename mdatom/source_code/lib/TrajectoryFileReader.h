#ifndef TRAJECTORYFILEREADER_H
#define TRAJECTORYFILEREADER_H

#include "MDParameters.h"

#include <string>
#include <vector>

class TrajectoryFileReader {
    
public:
    explicit TrajectoryFileReader(const std::string& inputFileName, const MDParameters& par);
    int getCalculationType() const { return calculationType; }
    std::size_t getNTimeframes() const { return nTimeframes; }
    const std::vector<std::vector<double>>& getVelocities() const { return velocities; }
private:
    int calculationType;
    std::size_t nAtoms;
    std::size_t nTimeframes;
    std::string fileName;
    std::string Title;
    std::vector<std::vector<double>> velocities;
};

#endif // TRACETORYFILEREADER_H
