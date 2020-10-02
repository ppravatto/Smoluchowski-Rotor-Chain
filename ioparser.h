#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <bits/stdc++.h> 
#include <fstream>
#include <string>
#include <sstream>
#include "math_utils.h"

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

namespace input{

    const int NUM_COMM = 18;                    //Total number of available commands
    std::string COMMAND_LIST[NUM_COMM] ={       //Commands list
        "NROT",
        "NUM-MIN",
        "ISOLATED-BASIS",
        "COUPLED-BASIS",
        "BARRIER",
        "DIFFUSION-COEFF",
        "NPT-INTEG-SINGLE",
        "ABS-INTEG-SINGLE",
        "REL-INTEG-SINGLE",
        "QAG-KEY-SINGLE",
        "NPT-INTEG-COUPLED",
        "ABS-INTEG-COUPLED",
        "REL-INTEG-COUPLED",
        "QAG-KEY-COUPLED",
        "SAVE-VQE",
        "SAVE-EIGVAL-LIST",
        "SINGLE-COUPLED-SCAN",
        "LOCKED-COUPLED-SCAN"
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
            bool init_flag, load_flag, vqe_key, eigval_list_key, scan_flag, locked_scan_flag;
            int num_rot, num_dihed, npt_int_single, npt_int_coupled, key_single, key_coupled;
            double abs_single, rel_single, abs_coupled, rel_coupled;
            int *basis_single, *basis_coupled, *num_mins, *scan_settings, *locked_scan_settings;
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
                npt_int_single = 10000; npt_int_coupled = 10000;
                key_single = 6; key_coupled = 6;
                abs_single = 1e-10; rel_single = 1e-10;
                abs_coupled = 1e-10; rel_coupled = 1e-10;
                vqe_key = false; eigval_list_key = false;
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
                scan_flag = false;
                locked_scan_flag = false;
            }

            ~INPUT_PARSER(){
                if(init_flag==true){
                    delete[] basis_single;
                    delete[] basis_coupled;
                    delete[] num_mins;
                    delete[] barrier;
                    delete[] diffusion;
                }
                if(scan_flag==true){
                    delete[] scan_settings;
                }
                if(locked_scan_flag==true){
                    delete[] locked_scan_settings;
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
                            case 6:
                                std::stringstream(line) >> npt_int_single;
                                break;
                            case 7:
                                std::stringstream(line) >> abs_single;
                                break;
                            case 8:
                                std::stringstream(line) >> rel_single;
                                break;
                            case 9:
                                std::stringstream(line) >> key_single;
                                break;
                            case 10:
                                std::stringstream(line) >> npt_int_coupled;
                                break;
                            case 11:
                                std::stringstream(line) >> abs_coupled;
                                break;
                            case 12:
                                std::stringstream(line) >> rel_coupled;
                                break;
                            case 13:
                                std::stringstream(line) >> key_coupled;
                                break;
                            case 14:
                                std::stringstream(line) >> vqe_key;
                                break;
                            case 15:
                                std::stringstream(line) >> eigval_list_key;
                                break;
                            case 16:
                                scan_flag = true;
                                scan_settings = new int[4];
                                read_list_line<int>(line, ",", scan_settings, 4);
                                if(scan_settings[3]==0){
                                    std::cout << "WARNING (input-parser): The scan step cannot be zero" << std::endl;
                                    std::cout << "                        -> Setting the step to 1" << std::endl;
                                    scan_settings[3] = 1;
                                }
                                break;
                            case 17:
                                locked_scan_flag = true;
                                locked_scan_settings = new int[3];
                                read_list_line<int>(line, ",", locked_scan_settings, 3);
                                if(locked_scan_settings[2]==0){
                                    std::cout << "WARNING (input-parser): The scan step cannot be zero" << std::endl;
                                    std::cout << "                        -> Setting the step to 1" << std::endl;
                                    locked_scan_settings[2] = 1;
                                }
                                break;
                            default:
                                break;
                        }
                    }
                }
                input.close();
                load_flag = true;
                if(scan_flag==true && locked_scan_flag==true){
                    std::cout << "WARNING (input-parser): Conflicting scan instructions were given in the input file" << std::endl;
                    std::cout << "                        -> Only the locked scan will be performed" << std::endl;
                    delete[] scan_settings;
                    scan_flag = false;
                }
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

            void copy_integrator_data(int* npt_int, int* key, double* abs, double* rel, bool single_flag){
                if(load_flag==false){
                    std::cout << "ERROR (input-parser): Input file not loaded yet" << std::endl;
                    exit(EXIT_FAILURE);
                }
                if(single_flag==true){
                    *npt_int = npt_int_single; *key=key_single;
                    *abs = abs_single; *rel = rel_single;
                }
                else{
                    *npt_int = npt_int_coupled; *key=key_coupled;
                    *abs = abs_coupled; *rel = rel_coupled;
                }
            }

            void get_general_settings(bool& vqe_key_, bool& eigval_list_key_, bool& scan_flag_, bool& locked_scan_flag_){
                vqe_key_ = vqe_key;
                eigval_list_key_ = eigval_list_key;
                scan_flag_ = scan_flag;
                locked_scan_flag_ = locked_scan_flag;
            }

            void coupled_scan_settings(int& dihedral_, int& start_, int& final_, int& step_){
                if(scan_flag==false){
                    std::cout << "ERROR (input-parser): Single dihedral scan has not been set" << std::endl;
                    exit(EXIT_FAILURE);
                }
                dihedral_ = scan_settings[0];
                start_ = scan_settings[1];
                final_ = scan_settings[2];
                step_ = scan_settings[3];
            }

            void coupled_locked_scan_settings(int& start_, int& final_, int& step_){
                if(locked_scan_flag==false){
                    std::cout << "ERROR (input-parser): Locked multiple dihedral scan has not been set" << std::endl;
                    exit(EXIT_FAILURE);
                }
                start_ = locked_scan_settings[0];
                final_ = locked_scan_settings[1];
                step_ = locked_scan_settings[2];
            }
    };
}

namespace output{

    std::string get_current_dir() {
        char buff[FILENAME_MAX];
        GetCurrentDir( buff, FILENAME_MAX );
        std::string current_working_dir(buff);
        return current_working_dir;
    }
    
}

#endif