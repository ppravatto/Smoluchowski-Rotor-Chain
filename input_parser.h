#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <bits/stdc++.h> 
#include <fstream>
#include <string>
#include <sstream>
#include "math_utils.h"

namespace input{

    const int NUM_COMM = 6;
    std::string COMMAND_LIST[NUM_COMM] ={
        "NROT",
        "NUM-MIN",
        "ISOLATED-BASIS",
        "COUPLED-BASIS",
        "BARRIER",
        "DIFFUSION-COEFF"
    };

    template <class T>
    void read_list_line(std::string line, std::string delimiter, T* vector, int length){
        int index = 0;
        size_t pos = 0;
        while ((pos = line.find(delimiter)) != std::string::npos) {
            std::string token = line.substr(0, pos);
            std::stringstream(token) >> vector[index];
            line.erase(0, pos + delimiter.length());
            index++;
            if(index>=length){
                std::cout << "ERROR (read_list_line): Too many arguments in line" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        std::stringstream(line) >> vector[index];
        if(index<length-1){
            std::cout << "ERROR (read_list_line): Too few arguments in line" << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    class INPUT_PARSER{
        private:
            std::string filename;
            bool init_flag, load_flag;
            int num_rot, num_dihed;
            int *basis_single, *basis_coupled, *num_mins;
            double *diffusion, *barrier;

        protected:
            bool check_command(std::string line){
                char first_char = line[0];
                return (first_char == '#')? true : false;
            }

            int get_command(std::string line){
                std::transform(line.begin(), line.end(), line.begin(), ::toupper);
                if(check_command(line)==false){
                    std::cout << "ERROR (input-parser): The line it is not marked as a command" << std::endl;
                    exit(EXIT_FAILURE);
                }
                int code = -1;
                for(int i=0; i<NUM_COMM; i++){
                    size_t pos = line.find(COMMAND_LIST[i]);
                    if(pos != std::string::npos){
                        code = i;
                        break;
                    }
                }
                return code;
            }

            void default_initializer(){
                if(init_flag==false){
                    std::cout << "ERROR (read_list_line): Memory buffers not yet allocated" << std::endl;
                    exit(EXIT_FAILURE);
                }
                for(int i=0; i<num_rot; i++){
                    diffusion[i] = 1.;
                    if(i<num_rot-1){
                        basis_single[i] = 100;
                        basis_coupled[i] = 10;
                        num_mins[i] = 1;
                        barrier[i] = 1.;
                    }
                }
            }

        public:
            INPUT_PARSER(std::string filename_){
                filename = filename_;
                init_flag = false;
                load_flag = false;
                std::ifstream input;
                input.open(filename);
                if(input.is_open()==false){
                    std::cout << """ERROR (input-parser): Failure in opening '""" + filename + """' file" << std::endl;
                    exit(EXIT_FAILURE);
                }

                while(!input.eof()){
                    std::string line;
                    std::getline(input, line);
                    if(check_command(line)==true){
                        if(get_command(line)==0){
                            std::getline(input, line);
                            std::stringstream(line) >> num_rot;
                            input.clear();
                            input.seekg(0);
                            break;
                        }
                    }
                }

                if(input.eof() || num_rot<=0){
                    std::cout << "ERROR (input-parser): Number of rotors not found or invalid" << std::endl;
                    exit(EXIT_FAILURE);
                }
                input.close();

                num_dihed = num_rot-1;
                try{
                    basis_single = new int [num_dihed];
                    basis_coupled = new int [num_dihed];
                    num_mins = new int [num_dihed];
                    barrier = new double [num_dihed];
                    diffusion = new double [num_rot];
                }
                catch(std::bad_alloc&){
                    std::cout << "ERROR (input-parser): Failure in allocating temporary memory buffer" << std::endl;
                    exit(EXIT_FAILURE);
                }
                init_flag = true;
            }

            ~INPUT_PARSER(){
                if(init_flag==true){
                    delete[] basis_single;
                    delete[] basis_coupled;
                    delete[] num_mins;
                    delete[] barrier;
                    delete[] diffusion;
                }
            }

            int load(){
                default_initializer();
                std::ifstream input;
                input.open(filename);
                while(!input.eof()){
                    std::string line;
                    std::getline(input, line);
                    if(check_command(line)==true){
                        int code = get_command(line);
                        std::getline(input, line);
                        switch (code){
                            case 1:
                                read_list_line<int>(line, ",", num_mins, num_dihed);
                                break;
                            case 2:
                                read_list_line<int>(line, ",", basis_single, num_dihed);
                                break;
                            case 3:
                                read_list_line<int>(line, ",", basis_coupled, num_dihed);
                                break;
                            case 4:
                                read_list_line<double>(line, ",", barrier, num_dihed);
                                break;
                            case 5:
                                read_list_line<double>(line, ",", diffusion, num_rot);
                                break;
                            default:
                                break;
                        }
                    }
                }
                input.close();
                load_flag = true;
                return num_rot;
            }

            void copy_system_data(int* num_mins_, int* basis_single_, int* basis_coupled_, double* barrier_, double* diffusion_){
                if(load_flag==false){
                    std::cout << "ERROR (input-parser): Input file not loaded yet" << std::endl;
                    exit(EXIT_FAILURE);
                }
                math_utils::copy_array<int>(num_mins, num_mins_, num_dihed);
                math_utils::copy_array<int>(basis_single, basis_single_, num_dihed);
                math_utils::copy_array<int>(basis_coupled, basis_coupled_, num_dihed);
                math_utils::copy_array<double>(barrier, barrier_, num_dihed);
                math_utils::copy_array<double>(diffusion, diffusion_, num_rot);
            }
    };

}

#endif