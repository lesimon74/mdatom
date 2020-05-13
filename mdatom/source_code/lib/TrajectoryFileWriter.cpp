#include "TrajectoryFileWriter.h"
#include "BinaryIO.h"
#include "CoordinatesAndVelocitiesInitializer.h" // For MAXTITLE value
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <utility>

using namespace std;

TrajectoryFileWriter::TrajectoryFileWriter(const MDParameters &parameters,
                                           std::string finalCoordFilename,
                                           std::string trajFilename,
                                           std::string velTrajFilename)
  : par(parameters),
    finalCoordinatesFilename(std::move(finalCoordFilename)),
    trajectoryCoordinatesFilename(std::move(trajFilename)),
    trajectoryVelocitiesFilename(std::move(velTrajFilename)) {
}

void TrajectoryFileWriter::writeBeforeRun() {
    ofstream fout1; // trajectory output
    ofstream fout2; // velocity trajectory output
    if (par.trajectoryOutput) {
        if (par.trajectoryOutputFormat == TrajectoryFileFormat::binary || par.trajectoryOutputFormat == TrajectoryFileFormat::ascii || par.trajectoryOutputFormat == TrajectoryFileFormat::positionAndVelocityBinary || par.trajectoryOutputFormat == TrajectoryFileFormat::positionAndVelocityAscii)
        {
            fout1.open(trajectoryCoordinatesFilename, ios::out);
            if (fout1.bad()) {
                throw std::runtime_error("can't open " + trajectoryCoordinatesFilename);
            }
            fout1 << par.title << endl;
        }
        if (par.trajectoryOutputFormat == TrajectoryFileFormat::velocityBinary || par.trajectoryOutputFormat == TrajectoryFileFormat::velocityAscii || par.trajectoryOutputFormat == TrajectoryFileFormat::positionAndVelocityBinary || par.trajectoryOutputFormat == TrajectoryFileFormat::positionAndVelocityAscii)
        {
            fout2.open(trajectoryVelocitiesFilename, ios::out); // deletes old file and creates new one from scratch
            if (fout2.bad()) {
                throw std::runtime_error("can't open " + trajectoryVelocitiesFilename);
            }
            fout2 << par.title << endl;
        }
    }
}

void TrajectoryFileWriter::writeFinalCoordinates(const std::vector<double>& positions,
                                                 const std::vector<double>& velocities) {
    if (par.finalXVOutput == FinalCoordinateFileFormat::ascii) {
        writeFinalCoordinatesInAsciiForm(positions, velocities);
    }
    else {
        writeFinalCoordinatesInBinaryForm(positions, velocities);
    }
}

void TrajectoryFileWriter::writeFinalCoordinatesInBinaryForm(const std::vector<double>& positions,
                                                             const std::vector<double>& velocities) {
    ofstream fout2;
    fout2.open(finalCoordinatesFilename, ios::out | ios::binary);
    if (fout2.bad()) {
        throw std::runtime_error("can't open " + finalCoordinatesFilename);
    }
    fout2.write(par.title.c_str(), MAXTITLE);
    BinaryIO::write(fout2, positions);
    BinaryIO::write(fout2, velocities);
}

void TrajectoryFileWriter::writeFinalCoordinatesInAsciiForm(const std::vector<double>& positions,
                                                            const std::vector<double>& velocities) {
    ofstream fout2;
    fout2.open(finalCoordinatesFilename, ios::out);
    if (fout2.bad()) {
        throw std::runtime_error("can't open " + finalCoordinatesFilename);
    }
    fout2 << par.title << "\n" ;
    fout2 << par.numberAtoms << "\n" ;
    for (int j = 0; j < par.numberAtoms; j++) {
        fout2 << setw(6) << j;
        for (int m = 0; m < 3; m++) {
            fout2 << setw(15) << positions[3 * j + m];
        }
        for (int m = 0; m < 3; m++) {
            fout2 << setw(15) << velocities[3 * j + m];
        }
        fout2 << "\n";
    }
}

void TrajectoryFileWriter::writeOutTrajectoryStep(const std::vector<double>& positions,
                                                  const std::vector<double>& velocities) {
    if (par.trajectoryOutput) {
        if (par.trajectoryOutputFormat == TrajectoryFileFormat::binary) {
            writeOutTrajectoryStepInBinaryForm(positions);
        } else if (par.trajectoryOutputFormat == TrajectoryFileFormat::ascii) {
            writeOutTrajectoryStepInAsciiForm(positions);
        } else if (par.trajectoryOutputFormat == TrajectoryFileFormat::velocityBinary) {
            writeOutTrajectoryStepVelocityInBinaryForm(velocities);
        } else if (par.trajectoryOutputFormat == TrajectoryFileFormat::velocityAscii) {
            writeOutTrajectoryStepVelocityInAsciiForm(velocities);
        } else if (par.trajectoryOutputFormat == TrajectoryFileFormat::positionAndVelocityBinary) {
            writeOutTrajectoryStepPosAndVelInBinaryForm(positions, velocities);
        } else if (par.trajectoryOutputFormat == TrajectoryFileFormat::positionAndVelocityAscii) {
            writeOutTrajectoryStepPosAndVelInAsciiForm(positions, velocities);
        }
    }
}

void TrajectoryFileWriter::writeOutTrajectoryStepInBinaryForm(const std::vector<double>& positions) {
    ofstream fileBW;
    fileBW.open(trajectoryCoordinatesFilename, ios::out | ios::app | ios::binary);
    if (fileBW.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryCoordinatesFilename);
    }
    BinaryIO::write(fileBW, positions);
}

void TrajectoryFileWriter::writeOutTrajectoryStepInAsciiForm(const std::vector<double>& positions) {
    ofstream fileFW;
    fileFW.open(trajectoryCoordinatesFilename, ios::out | ios::app);
    if (fileFW.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryCoordinatesFilename);
    }
    int ntot = 3 * par.numberAtoms;
    for (int i = 0; i < ntot; i += 10) {
        for (int j = i; (j < i + 10 && j < ntot); j++) {
            fileFW << setw(10) << positions[j];
        }
        fileFW << endl;
    }
}


void TrajectoryFileWriter::writeOutTrajectoryStepVelocityInBinaryForm(const std::vector<double>& velocities) {
    ofstream fileBW;
    fileBW.open(trajectoryVelocitiesFilename, ios::out | ios::app | ios::binary);
    if (fileBW.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryVelocitiesFilename);
    }
    BinaryIO::write(fileBW, velocities);
}


void TrajectoryFileWriter::writeOutTrajectoryStepVelocityInAsciiForm(const std::vector<double>& velocities) {
    ofstream fileFW;
    fileFW.open(trajectoryVelocitiesFilename, ios::out | ios::app);
    if (fileFW.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryVelocitiesFilename);
    }
    int ntot = 3 * par.numberAtoms;
    for (int i = 0; i < ntot; i += 10) {
        for (int j = i; (j < i + 10 && j < ntot); j++) {
            fileFW << setw(10) << velocities[j];
        }
        fileFW << endl;
    }
}


void TrajectoryFileWriter::writeOutTrajectoryStepPosAndVelInBinaryForm(const std::vector<double>& positions,
                                                                       const std::vector<double>& velocities) {
    ofstream fileBW;
    fileBW.open(trajectoryCoordinatesFilename, ios::out | ios::app | ios::binary);
    if (fileBW.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryCoordinatesFilename);
    }
    BinaryIO::write(fileBW, positions);
    
    ofstream fileBW2;
    fileBW2.open(trajectoryVelocitiesFilename, ios::out | ios::app | ios::binary);
    if (fileBW2.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryVelocitiesFilename);
    }
    BinaryIO::write(fileBW2, velocities);
}


void TrajectoryFileWriter::writeOutTrajectoryStepPosAndVelInAsciiForm(const std::vector<double>& positions,
                                                                      const std::vector<double>& velocities) {
    ofstream fileFW;
    fileFW.open(trajectoryCoordinatesFilename, ios::out | ios::app);
    if (fileFW.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryCoordinatesFilename);
    }
    int ntot = 3 * par.numberAtoms;
    for (int i = 0; i < ntot; i += 10) {
        for (int j = i; (j < i + 10 && j < ntot); j++) {
            fileFW << setw(10) << positions[j];
        }
        fileFW << endl;
    }
    
    ofstream fileFW2;
    fileFW2.open(trajectoryVelocitiesFilename, ios::out | ios::app);
    if (fileFW2.bad()) {
        throw runtime_error("I/O ERROR: cannot write to file: " + trajectoryVelocitiesFilename);
    }
    int ntot2 = 3 * par.numberAtoms;
    for (int i = 0; i < ntot2; i += 10) {
        for (int j = i; (j < i + 10 && j < ntot2); j++) {
            fileFW2 << setw(10) << velocities[j];
        }
        fileFW2 << endl;
    }
}
