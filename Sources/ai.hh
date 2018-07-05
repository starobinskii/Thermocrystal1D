#pragma once

#include <array>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <sys/stat.h>

namespace ai{
    inline std::string getVersion(){
        return "1.1.0";
    }

    template<typename T>
    std::string string(const T value){
        std::ostringstream stream;

        stream << value;

        return stream.str();
    }

    inline bool hasPrefix(const std::string &str, const std::string &prefix){
        return str.size() >= prefix.size()
            && 0 == str.compare(0, prefix.size(), prefix);
    }

    inline bool hasSuffix(const std::string &str, const std::string &suffix){
        return str.size() >= suffix.size() &&
            str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    inline bool assignStringParameter(
        const char *input,
        const std::string name,
        std::string &value
    ){
        std::string parameter = std::string(input);

        if(ai::hasPrefix(parameter, name)){
            value = parameter.substr(name.size());

            return true;
        }

        return false;
    }

    inline bool assignBooleanParameter(
        const char *input,
        const std::string name,
        bool &value
   ){
       if(name == std::string(input)){
           value = true;

           return true;
       }

       return false;
   }

    template<typename T>
    inline bool assignParameter(
        const char *input,
        const std::string name,
        T &value
    ){
        std::string parameter = std::string(input);

        if(ai::hasPrefix(parameter, name)){
            parameter = parameter.substr(name.size());
            
            if(std::istringstream(parameter) >> value){
                return true;
            }
        }

        return false;
    }

    inline void showProgressBar(double progress){
        if(1 < progress){
            progress = 1;
        }

        if(0.01 > progress){
            progress = 0;
        }

        int width = progress * 73;

        std::cout << std::fixed;
        std::cout.precision(1);
        std::cout.flush();

        std::cout << "\r" << std::string(width, '=')
            << std::string(73 - width, '-')
            << " " << progress * 100. << "%";
    }

    template<typename T>
    inline void saveMatrix(
        const std::string filename,
        const std::vector<std::vector <T> > matrix,
        std::string comment = std::string(),
        std::string type = std::string("text"),
        std::string delimiter = std::string(" ")
    ){
        std::string extension("_v.txt");
        std::string prefix("");
        std::string suffix("");

        if(std::string("wolfram") == type){
            extension = std::string("_v.wm");
            prefix = std::string("{");
            delimiter = std::string(", ");
            suffix = std::string("}");
        }

        if(std::string("excel") == type){
            extension = std::string("_v.csv");
            delimiter = std::string("; ");
        }

        if(std::string("data") == type){
            extension = std::string("_v.dat");
            delimiter = std::string("\t");
        }
        
        std::ofstream output(filename + extension);

        if(!output.good()){
            throw std::runtime_error(
                ai::string("Exception while saving the matrix into the file: ") 
                + filename
            );
        }

        if(std::string() != comment){
            output << comment << std::endl;
        }

        output << prefix;

        for(const std::vector<T> &row : matrix){
            output << prefix;
            
            const std::size_t lastIndex = row.size() - 1;

            for(std::size_t i = 0; i < lastIndex; ++i){
                output << std::setw(14) << row[i] << delimiter;
            }

            output << std::setw(14) << row[lastIndex] << suffix << std::endl;
        }
        
        output << suffix;
    }

    template<typename T>
    inline void saveVector(
        const std::string filename,
        const std::vector <T> vector,
        std::string comment = std::string(),
        std::string type = std::string("text"),
        std::string delimiter = std::string("\n")
    ){
        std::string extension("_v.txt");
        std::string prefix("");
        std::string suffix("");

        if(std::string("wolfram") == type){
            extension = std::string("_v.wm");
            prefix = std::string("{");
            delimiter = std::string(", ");
            suffix = std::string("}");
        }

        if(std::string("excel") == type){
            extension = std::string("_v.csv");
            delimiter = std::string("; ");
        }

        if(std::string("data") == type){
            extension = std::string("_v.dat");
            delimiter = std::string("\t");
        }

        std::ofstream output(filename + extension);

        if(!output.good()){
            throw std::runtime_error(
                ai::string("Exception while saving the vector into the file: ") 
                + filename
            );
        }

        if(std::string() != comment){
            output << comment << std::endl;
        }

        output << prefix;

        const std::size_t lastIndex = vector.size() - 1;

        for(std::size_t i = 0; i < lastIndex; ++i){
            output << vector[i] << delimiter;
        }

        output << vector[lastIndex] << suffix;
    }
}