#include "CorrelationCalculator.h"

#include <fftw3.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

CorrelationCalculator::CorrelationCalculator(std::unique_ptr<TrajectoryFileReader> reader): nTimeframes{reader->getNTimeframes()}, velocities{reader->getVelocities()}, fReader{std::move(reader)}
{
    nAtoms = velocities.at(0).size();
    correlations.resize(nTimeframes);

    // Initialize meanVelocities
    meanVelocities.resize(nAtoms);
    for (int i = 0; i < nAtoms; ++i) {
        double sum{0};
        for (int j = 0; j < nTimeframes; ++j) {
            sum += velocities.at(j).at(i);
        }
        sum /= nTimeframes;
        meanVelocities.at(i) = sum;
    }

    // Initialize std. Dev.
    standardDeviationPerAtom.resize(nAtoms);
    for (int i = 0; i < nAtoms; ++i) {
        double sum{0};
        for (int j = 0; j < nTimeframes; ++j) {
            sum += std::pow(velocities.at(j).at(i) - meanVelocities.at(i), 2);
        }
        sum /= nTimeframes;
        standardDeviationPerAtom.at(i) = std::sqrt(sum);
    }
}

void CorrelationCalculator::printCorrelation()
{
    calculateCorrelation();
    for (std::size_t i = 0; i < correlations.size(); ++i) {
        std::cout << i << std::setw(15) << correlations.at(i) << std::endl;
    }
}

void CorrelationCalculator::calculateCorrelation()
{
    auto start_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();
    
    if (fReader->getCalculationType() == 0) { // direct multiplication method
        auto start_time = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < nTimeframes; ++i) {
            for (std::size_t j = 0; j < nAtoms; ++j) {
                double sum{0};
                for (std::size_t k = 0; k < nTimeframes-i; ++k) {
                    sum += (velocities.at(k).at(j)-meanVelocities.at(j))*(velocities.at(k+i).at(j)-meanVelocities.at(j));
                }
                correlations[i] += (sum/(nTimeframes-i)/standardDeviationPerAtom.at(j)/standardDeviationPerAtom.at(j));
            }
            correlations[i] /= nAtoms;
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        auto time = end_time - start_time;
        std::cout << "Time required to calculate correlation: " << time/std::chrono::milliseconds(1) << " ms.\n";
    } else if (fReader->getCalculationType() == 1) { // FFT method
        auto start_time = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < nAtoms; ++i) {
            fftw_complex *in, *out;
            int N = static_cast<int>(nTimeframes)*2; // due to zero filling.
            in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
            out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
            for (int j = 0; j < nTimeframes; ++j) {
                in[j][0] = velocities.at(j).at(i)-meanVelocities.at(i);
                in[j+nTimeframes][0] = in[j][1] = in[j+nTimeframes][1] = out[j][0] = out[j][1] = out[j+nTimeframes][0] = out[j+nTimeframes][1] = 0;
            }
            fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

            fftw_execute(p);
            for (int j = 0; j < 2*nTimeframes; ++j) {
                // multiplication with complex conjugate
                in[j][0] = out[j][0]*out[j][0] + out[j][1]*out[j][1];
                in[j][1] = out[j][0] = out[j][1] = 0;
            }
            // contrary to what the name suggests, this is a back transform.
            fftw_plan p2 = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p2);
            
            for (int j = 0; j < nTimeframes; ++j) {
                // Comment out the 2 lines below to use absolute value of FT.
                // double norm = std::sqrt(out[j][0]*out[j][0] + out[j][1]*out[j][1]);
                // correlations[j] += (norm/(2*nTimeframes*(nTimeframes-j)));
                correlations[j] += ((out[j][0]/(2*nTimeframes*(nTimeframes-j)))/standardDeviationPerAtom.at(i)/standardDeviationPerAtom.at(i));
            }
            fftw_free((fftw_complex*)in);
            fftw_free((fftw_complex*)out);
        }
        for (int j = 0; j < nTimeframes; ++j) {
            correlations[j] /= nAtoms;
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        auto time = end_time - start_time;
        std::cout << "Time required to calculate correlation: " << time/std::chrono::milliseconds(1) << " ms.\n";
    } else if (fReader->getCalculationType() == 2) { // Calculate correlation only for a specific atom.
        // Use std::cerr so it doesn't get redirected to file if user tries to do so. (e.g. ./mdatom params.inp > out.txt)
        std::cerr << "Please enter the number of the atom for which you would like to calculate the correlation (1-" << nAtoms << "): " << std::flush;
        // Assume the user is not retarded and a valid input is given.
        double i{0};
        std::cin >> i;
        fftw_complex *in, *out;
        int N = static_cast<int>(nTimeframes)*2; // due to zero filling.
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
        for (int j = 0; j < nTimeframes; ++j) {
            in[j][0] = velocities.at(j).at(i)-meanVelocities.at(i);
            in[j+nTimeframes][0] = in[j][1] = in[j+nTimeframes][1] = out[j][0] = out[j][1] = out[j+nTimeframes][0] = out[j+nTimeframes][1] = 0;
        }
        fftw_plan p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

        fftw_execute(p);
        for (int j = 0; j < 2*nTimeframes; ++j) {
            // multiplication with complex conjugate
            in[j][0] = out[j][0]*out[j][0] + out[j][1]*out[j][1];
            in[j][1] = out[j][0] = out[j][1] = 0;
        }
        // contrary to what the name suggests, this is a back transform.
        fftw_plan p2 = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(p2);

        for (int j = 0; j < nTimeframes; ++j) {
            // Comment out the 2 lines below to use absolute value of FT.
            // double norm = std::sqrt(out[j][0]*out[j][0] + out[j][1]*out[j][1]);
            // correlations[j] += (norm/(2*nTimeframes*(nTimeframes-j)));
            correlations[j] += ((out[j][0]/(2*nTimeframes*(nTimeframes-j)))/standardDeviationPerAtom.at(i)/standardDeviationPerAtom.at(i));
        }
        fftw_free((fftw_complex*)in);
        fftw_free((fftw_complex*)out);
    } else {
        std::cout << "Please set XVOutput to 0 (direct multiplication method), 1 (FFT method) or 2 (for calculating the correlation only of a specific atom)\n";
    }
}
