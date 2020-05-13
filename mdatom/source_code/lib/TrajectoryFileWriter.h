#ifndef TRAJECTORYFILEWRITER_H
#define TRAJECTORYFILEWRITER_H

#include "MDParameters.h"
#include <string>
#include <vector>

/*!
 * This class is responsible for writing out the atom coordinates during a MD simulation.
 */
class TrajectoryFileWriter {
  public:
    TrajectoryFileWriter(const MDParameters& parameters,
                         std::string finalCoordFilename,
                         std::string trajFilename,
                         std::string velTrajFilename);
    void writeBeforeRun();
    void writeFinalCoordinates(const std::vector<double>& positions, const std::vector<double>& velocities);
    void writeOutTrajectoryStep(const std::vector<double>& positions, const std::vector<double>& velocities);

  private:
    void writeFinalCoordinatesInBinaryForm(const std::vector<double>& positions, const std::vector<double>& velocities);
    void writeFinalCoordinatesInAsciiForm(const std::vector<double>& positions, const std::vector<double>& velocities);
    void writeOutTrajectoryStepInBinaryForm(const std::vector<double>& positions);
    void writeOutTrajectoryStepInAsciiForm(const std::vector<double>& positions);
    void writeOutTrajectoryStepVelocityInBinaryForm(const std::vector<double>& velocities);
    void writeOutTrajectoryStepVelocityInAsciiForm(const std::vector<double>& velocities);
    void writeOutTrajectoryStepPosAndVelInBinaryForm(const std::vector<double>& positions,
                                                     const std::vector<double>& velocities);
    void writeOutTrajectoryStepPosAndVelInAsciiForm(const std::vector<double>& positions,
                                                    const std::vector<double>& velocities);

    const MDParameters& par;
    const std::string finalCoordinatesFilename;
    const std::string trajectoryCoordinatesFilename;
    const std::string trajectoryVelocitiesFilename;
};

#endif // TRAJECTORYFILEWRITER_H
