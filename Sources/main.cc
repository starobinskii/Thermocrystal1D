#include <iostream>
#include "mpi.h"
#include "ai.hh"

#include "thermocrystal1D.hh"

int main(int argc, char *argv[]){
    bool standingWave = false;
    bool runningWave = false;
    bool saveVelocities = false;
    
    size_t id = 0;
    size_t length = 1000;
    
    double dispersion = 1.;
    
    std::string path("./");    
    std::string filename("");
    
    for(int i = 1; i < argc; ++i){
        if("-h" == std::string(argv[i]) || "--help" == std::string(argv[i])){
            std::cout << "usage: averageMatrix [options]"
                << std::endl
                << "    -h  --help            print this usage and exit"
                << std::endl<< std::endl
                
                << "    --path=<path>         path to save files [string]"
                << std::endl
                << "    --filename=<path>     name of the output file [string]"
                << std::endl
                << "    --id=<value>          positive identification number "
                << "[size_t]"
                << std::endl
                << "    --length=<value>      length of the crystal in "
                << "particles [size_t]"
                << std::endl
                << "    --dispersion=<value>  ratio of thermal energy to "
                << "mechanic energy [double]"
                << std::endl
                << "    --save-velocities     flag to save velocities [bool]"
                << std::endl
                << "    --standing-wave       flag to calculate a stanging "
                << "wave [bool]"
                << std::endl
                << "    --running-wave        flag to calculate a running "
                << "wave [bool]"
                << std::endl;
            
            return 0;
        }
        
        if(
            ai::assignStringParameter(argv[i], "--path=", path)
            || ai::assignStringParameter(argv[i], "--filename=", filename)
            || ai::assignParameter(argv[i], "--id=", id)
            || ai::assignParameter(argv[i], "--length=", length)
            || ai::assignParameter(argv[i], "--dispersion=", dispersion)
            || ai::assignBooleanParameter(argv[i], "--save-velocities", 
                saveVelocities)
            || ai::assignBooleanParameter(argv[i], "--standing-wave", 
                standingWave)
            || ai::assignBooleanParameter(argv[i], "--running-wave", 
                runningWave)
        ){
            continue;
        }
    }
    
    if(!standingWave && !runningWave){
        std::cout << "No wave has been specified. Go with a running wave."
            << std::endl;
        runningWave = true;
    }
    
    if(!ai::hasSuffix(path, "/")){
        path += std::string("/");
    }
    
    path += filename;
    
    if(0 < id){
        filename += ai::string(id);
    }
    
    MPI::Init();
    
    if(runningWave){
        thermocrystal1D(
            path + std::string("rw"),
            std::string("running"),
            dispersion,
            length, 
            saveVelocities
        );
    }
    
    MPI::COMM_WORLD.Barrier();
    if(standingWave){
        thermocrystal1D(
            path + std::string("sw"),
            std::string("standing"),
            dispersion,
            length,
            saveVelocities
        );
    }
    MPI::Finalize();
    
    return 0;
}
