#include "TrajectoryFileReader.h"

#include <iostream> // To remove later
#include <cmath>
#include <fstream>

TrajectoryFileReader::TrajectoryFileReader(const std::string& inputFileName, const MDParameters& par):
    calculationType{trajectoryFileFormatToInt(par.trajectoryOutputFormat)},
    nAtoms{static_cast<std::size_t>(par.numberAtoms)},
    nTimeframes{static_cast<std::size_t>(par.numberMDSteps/par.trajectoryOutputInterval)},
    fileName{inputFileName}
{
    std::ifstream fin;
    fin.open(fileName, std::ios::in);
        if (fin.bad()) {
            throw std::runtime_error("can't open " + fileName);
        }
    getline(fin,Title);

    double x, y, z;
    // Assuming a velocity trajectory file was read in, not a coordinate trajectory file.
    velocities.resize(nTimeframes);
    for (std::size_t i = 0; i < nTimeframes; ++i) {
        velocities.at(i).resize(nAtoms);
        // Read in the velocities
        for (std::size_t j = 0; j < nAtoms; ++j) {
            fin >> x;
            fin >> y;
            fin >> z;
            // Save absolute value of velocity instead of its x, y and z components.
            velocities.at(i).at(j) = std::sqrt(std::pow(x, 2) + std::pow(y, 2) + std::pow(z, 2));
        }
    }
}
