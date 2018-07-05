#pragma once

void applyPeriodic(std::vector<double> &quantity);

/*/
 * dispersion is a ratio of thermal energy to mechanic energy
 * set waveType either to "running" or to "standing"
/*/
void thermocrystal1D(
    const std::string filename,
    const std::string waveType,
    const double dispersion = 1.,
    const double maximumTime = 100000,
    size_t length = 1000,
    const bool saveVelocities = false,
    const double dt = 0.05,
    const bool runningFromCLI = false
);