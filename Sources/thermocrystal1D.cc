#include <iostream>
#include <fstream>
#include <cstdlib>
#include <random>
//#include "mpi.h"
#include "ai.hh"

#if !defined(M_PI)
    #define M_PI 3.141593
#endif

void applyPeriodic(std::vector<double> &quantity){
    const size_t length = quantity.size();
    
    quantity[0] = quantity[length - 2];
    quantity[length - 1] = quantity[1];
}

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
){
    
    bool runningWaveFlag = false;
    
    const size_t stepToSaveData = 3000;
    const size_t precision = 100000000;
    //const size_t rank = (size_t) MPI_COMM_WORLD.Get_rank();
    //const size_t size = (size_t) MPI_COMM_WORLD.Get_size();
    
    const double dx = 1.;
    const double alpha = 1.;
    double amplitude = 0.02;
    const double omegaSquare = 1.;
    
    if(!runningFromCLI){
        std::cout << "Incoming parameters:" << std::endl
            << waveType << " wave with " << length << " particles," << std::endl
            << "dispersion is " << dispersion << ", time is " << maximumTime 
            << std::endl << "cluster size is " << 1 << std::endl;
    }
    
    std::vector<double> displacement(length + 2, 0.);
    std::vector<double> velocity(length + 2, 0.);
    
    if(std::string("running") == waveType){
        runningWaveFlag = true;
        /*/ running wave goes both way /*/
        /*/ we calculate propagation only in one direction /*/
        amplitude = 0.01;
    }
    
    size_t step = 0;
    
    double currentTime = 0.;
    /*/ need in case of a running wave only /*/
    const double shiftTime = dx / sqrt(omegaSquare);
    double savedTime = 0.;
    
    // double *arrayBuffer;
    // double *velocityArray;
    
    // double valueBuffer = 0;
    
    std::vector<double> energyOutput;
    
    std::vector< std::vector<double> > velocityOutput;
    
    // if(saveVelocities){
    //     velocityArray = new double[length];
        
    //     if(0 == rank){
    //         arrayBuffer = new double[size * length];
    //     }
    // }
    
    std::random_device rd;
    srand(rd() * time(0));
    
    /*/ now applying initial conditions /*/
    double sum = 0;
    
    /*/ thermal velocities /*/
    for(size_t i = 1; i <= length; ++i){
        velocity[i] = amplitude * sqrt(dispersion * 3. / 2.) 
            * (2. / precision) * ((rand() % precision) - 1.);
        sum += velocity[i];
    }
    
    sum /= (double) length;
    
    /*/ mechanic wave velocities /*/
    for(size_t i = 1; i <= length; ++i){
        double x = 2. * M_PI * (i - 1.) / ((double) length);
        
        velocity[i] += amplitude * sin(x) - sum;
        
        /*/ for a case of a standing wave displacements are already zero /*/
        if(runningWaveFlag){
            displacement[i] = amplitude * cos(x) 
                * ((double) length) / (2. * M_PI * sqrt(omegaSquare));
        }
    }
    
    /*/ now applying periodic boundary conditions /*/
    applyPeriodic(velocity);
    applyPeriodic(displacement);
    
    while(maximumTime > currentTime){
        /*/ calculating a step /*/
        for(size_t i = 1; i <= length; ++i){
            velocity[i] += (displacement[i + 1] - 2 * displacement[i] 
                + displacement[i - 1]) * (omegaSquare + alpha * (displacement[i + 1] 
                - displacement[i - 1])) * dt / dx / dx;
        }
        
        for(size_t i = 1; i <= length; ++i){
            displacement[i] += velocity[i] * dt;
        }
        
        applyPeriodic(velocity);
        applyPeriodic(displacement);
        
        currentTime += dt;
        ++step;
        
        /*/ save the data /*/
        if(stepToSaveData < step){
            step = 0;
            
            // if(saveVelocities){
            //     /*/ gather the velocities /*/
            //     for(size_t i = 1; i <= length; ++i){
            //         velocityArray[i - 1] = velocity[i];
            //     }
                
            //     MPI_COMM_WORLD.Reduce(velocityArray, arrayBuffer, 
            //         (int) length, MPI_DOUBLE, MPI_SUM, 0);
                
            //     if(0 == rank){
            //         std::vector<double> velocityBuffer;
                    
            //         for(size_t i = 0; i < length; ++i){
            //             velocityBuffer.push_back(arrayBuffer[i]);
            //         }
                    
            //         velocityOutput.push_back(velocityBuffer);
            //     }
            // }
            
            /*/ calculating energy /*/
            double kineticEnergy = 0;
            double potentialEnergy = 0;
            
            for(int i = 1; i <= length; ++i){
                float x = (i - 1.) / ((double) length);
                
                kineticEnergy += velocity[i] * sin(2 * M_PI * x);
                potentialEnergy += (displacement[i - 1] - displacement[i]) 
                    * cos(2 * M_PI *x);
            }
            kineticEnergy *= kineticEnergy;
            potentialEnergy *= potentialEnergy;
            
            double energy = potentialEnergy + kineticEnergy;
            
            /*/ gather the energy /*/
           // MPI_COMM_WORLD.Reduce(&energy, &valueBuffer, 1, MPI_DOUBLE, MPI_SUM, 0);
            
            //if(0 == rank){
                energyOutput.push_back(energy);
                
                if(!runningFromCLI){
                    ai::showProgressBar(currentTime / maximumTime);
                }
            //}
        }
        
        /*/ in cave of a running wave move the wave forward /*/
        if(runningWaveFlag && shiftTime < currentTime - savedTime){
            for(size_t i = 1; i <= length; ++i){
                displacement[i] = displacement[i + 1];
                velocity[i] = velocity[i + 1];
            }
            
            applyPeriodic(velocity);
            applyPeriodic(displacement);
            
            savedTime += shiftTime;
        }
    }
    
    if(!runningFromCLI){
        ai::showProgressBar(1.);
    }
    
    /*/ saving the data /*/
    // if(0 == rank){
    //     if(saveVelocities){
    //         ai::saveMatrix(filename, velocityOutput);
    //     }
        
    ai::saveVector(filename, energyOutput);
    // }
    
    // if(saveVelocities){
    //     delete [] velocityArray;
        
    //     if(0 == rank){
    //         delete [] arrayBuffer;
    //     }
    // }
}
