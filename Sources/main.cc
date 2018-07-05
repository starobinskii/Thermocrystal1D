#include <iostream>
#include "mpi.h"
#include "ai.hh"

#include "thermocrystal1D.hh"

int main(int argc, char *argv[]){
    bool standingWave = false;
    bool runningWave = false;
    bool saveVelocities = false;
    bool runningFromCLI = false;
    
    size_t id = 0;
    size_t length = 1000;
    
    double dispersion = 1.;
    double modelingTime = 100000;
    
    std::string path("./");    
    std::string filename("");
    
    for(int i = 1; i < argc; ++i){
        if("-h" == std::string(argv[i]) || "--help" == std::string(argv[i])){
            std::cout << "usage: ./thermocrystal1D [options]"
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
                << "    --time=<value>        modeling time [double]"
                << std::endl
                << "    --save-velocities     flag to save velocities [bool]"
                << std::endl
                << "    --standing-wave       flag to calculate a stanging "
                << "wave [bool]"
                << std::endl
                << "    --running-wave        flag to calculate a running "
                << "wave [bool]"
                << std::endl
                << "    --env=cli             flag for CLI execution [bool]"
                << std::endl;
            
            return 0;
        }
        
        if(
            ai::assignStringParameter(argv[i], "--path=", path)
            || ai::assignStringParameter(argv[i], "--filename=", filename)
            || ai::assignParameter(argv[i], "--id=", id)
            || ai::assignParameter(argv[i], "--length=", length)
            || ai::assignParameter(argv[i], "--dispersion=", dispersion)
            || ai::assignParameter(argv[i], "--time=", modelingTime)
            || ai::assignBooleanParameter(argv[i], "--save-velocities", 
                saveVelocities)
            || ai::assignBooleanParameter(argv[i], "--standing-wave", 
                standingWave)
            || ai::assignBooleanParameter(argv[i], "--running-wave", 
                runningWave)
            || ai::assignBooleanParameter(argv[i], "--env=cli", 
                runningFromCLI)
        ){
            continue;
        }
    }
    
    if(!standingWave && !runningWave){
        if(!runningFromCLI){
            std::cout << "No wave has been specified. Go with a running wave."
                << std::endl;
        }
        runningWave = true;
    }
    
    if(!ai::hasSuffix(path, "/")){
        path += std::string("/");
    }
    
    const std::string folder(path);
    
    path += filename;
    
    MPI::Init();
    
    const size_t rank = (size_t) MPI::COMM_WORLD.Get_rank();
    const size_t size = (size_t) MPI::COMM_WORLD.Get_size();
    
    if(runningWave){
        thermocrystal1D(
            path + ai::string("rw_") + ai::string(id),
            std::string("running"),
            dispersion,
            modelingTime,
            length,
            saveVelocities,
            0.05,
            runningFromCLI
        );
    }
    
    MPI::COMM_WORLD.Barrier();
    
    if(standingWave){
        thermocrystal1D(
            path + ai::string("sw") + ai::string(id),
            std::string("standing"),
            dispersion,
            modelingTime,
            length,
            saveVelocities,
            0.05,
            runningFromCLI
        );
    }
    
    if(0 == rank){
        ai::saveVector(
            folder + ai::string("output") + ai::string(id),
            std::vector<double>{(double) length, dispersion, (double) size},
            std::string("#Length; dispersion; cluster size")
        );
    }
    
    MPI::Finalize();
    
    return 0;
}
