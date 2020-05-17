#include "CorrelationCalculator.h"

#include <fftw3.h>

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
    if (fReader->getCalculationType() == 0) { // direct multiplication method
        for (std::size_t i = 0; i < nTimeframes; ++i) {
            for (std::size_t j = 0; j < nAtoms; ++j) {
                double sum{0};
                for (std::size_t k = 0; k < nTimeframes-i; ++k) {
                    sum += velocities.at(i).at(j)*velocities.at(k+i).at(j);
                }
                correlations[i] += (sum/(nTimeframes-i));
            }
            correlations[i] /= nAtoms;
        }
    } else if (fReader->getCalculationType() == 1) { // FFT method
        for (std::size_t i = 0; i < nAtoms; ++i) {
            fftw_complex *in, *out;
            int N = static_cast<int>(nTimeframes)*2; // due to zero filling.
            in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
            out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
            for (int j = 0; j < nTimeframes; ++j) {
                in[j][0] = velocities.at(j).at(i);
                in[j+nTimeframes][0] = in[j][1] = in[j+nTimeframes][1] = 0;
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
                correlations[j] += (out[j][0]/(2*nTimeframes*(nTimeframes-j)));
            }
        }
        for (int j = 0; j < nTimeframes; ++j) {
            correlations[j] /= nAtoms;
        }
    } else {
        std::cout << "Please set XVOutput to 0 (direct multiplication method) or 1 (FFT method)\n";
    }
}
