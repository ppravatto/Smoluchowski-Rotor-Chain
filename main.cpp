#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string>
#include <limits.h>
#include <omp.h>
#include "rotors.h"
#include "ioparser.h"

using namespace std;

int main(int argc, char** argv){

    if(argc==1){
        cout << "ERROR (main.cpp): No input file selected" << endl;
        exit(EXIT_FAILURE);
    }
    string filename = argv[1];
    cout << "Selected input file: " << filename << endl << endl;

    if(argc>2){
        cout << "WARNING (main.cpp): Too many arguments received, only the first will be considered" << endl;
    }

    string base_dir = output::get_current_dir();

    input::INPUT_PARSER datafile(filename);             //Load input file with the input parser
    int N_rot = datafile.load();                        //Set total number of rotors
    int N = N_rot-1;                                    //Number od dihedral angles

    bool vqe_key = false, eigval_save_key = false, scan_flag = false, locked_scan_flag = false, debug_key = false;
    datafile.get_general_settings(vqe_key, eigval_save_key, scan_flag, locked_scan_flag, debug_key);

    if(debug_key == true){
        cout << "                          *** DEBUG MODE ***" << endl;
    }

    double * D = new double [N_rot];                    //Allocate an array to store the diffusion coefficients for each rotor
    double * dihedral_barrier = new double [N];         //Allocate an array to store the barrier height for each dihedral
    int * dihedral_num_mins = new int [N];              //Allocate an array to store the type of potential function for each dihedral
    int * single_dihedral_basis_order = new int [N];    //Allocate an array to store the number of basis functions to use for generating each dihedral basis set
    int * composite_basis_set_cutoff = new int [N];     //Allocate an array to store the total number of single rotor basis functions to be use for generating the composit basis set

    //Retrieve system data from input parser 
    datafile.copy_system_data(
                                dihedral_num_mins,
                                single_dihedral_basis_order,
                                composite_basis_set_cutoff,
                                dihedral_barrier,
                                D
                            );

    int scan_index=-1, start_scan=0, end_scan=0, step_scan=0;
    if(scan_flag==true){
        datafile.coupled_scan_settings(scan_index, start_scan, end_scan, step_scan);
        if(start_scan > single_dihedral_basis_order[scan_index] || end_scan > single_dihedral_basis_order[scan_index]){
            cout << "ERROR (main.cpp): The selected single dihedral basis set is smaller than the required scan range" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(locked_scan_flag==true){
        datafile.coupled_locked_scan_settings(start_scan, end_scan, step_scan);
        int min_single_basis = single_dihedral_basis_order[0];
        for(int i=0; i<N; i++){
            if(single_dihedral_basis_order[i] < min_single_basis){
                min_single_basis = single_dihedral_basis_order[i];
            }
        }
        if(start_scan > min_single_basis || end_scan > min_single_basis){
            cout << "ERROR (main.cpp): The selected single dihedral basis set is smaller than the required scan range" << endl;
            exit(EXIT_FAILURE);
        }
    }

    if((scan_flag==true || locked_scan_flag==true) && start_scan==end_scan){
        cout << "WARNING (main.cpp): The start point in a scan instruction cannot be equal to the end point" << endl;
        cout << "                    -> Falling back to regular single point calculation" << endl;
        scan_flag = false;
        locked_scan_flag = false;
    }

    //Retrieve integration parameters for the single dihedral from input parser 
    int npt_int_single, key_int_single;
    double abs_int_single, rel_int_single;
    datafile.copy_integrator_data(&npt_int_single, &key_int_single, &abs_int_single, &rel_int_single, true);

    //Retrieve integration parameters for the single dihedral from input parser 
    int npt_int_system, key_int_system;
    double abs_int_system, rel_int_system;
    datafile.copy_integrator_data(&npt_int_system, &key_int_system, &abs_int_system, &rel_int_system, false);

    for(int i=0; i<N; i++){
        if(composite_basis_set_cutoff[i]>single_dihedral_basis_order[i]){
            cout << "ERROR (main.cpp): The cutoff for the dihedral " << i << " cannot be smaller than the original basis set" << endl;
            exit(EXIT_FAILURE); 
        }
    }

    int M = 1;                                          //Total number of basis functions
    for(int i=0; i < N; i++){
        if(composite_basis_set_cutoff[i] >= int(INT_MAX/M)){
            cout << "ERROR (main.cpp): Integer overflow in generating the total number of basis functions" << endl;
            exit(EXIT_FAILURE);
        }
        M *= composite_basis_set_cutoff[i];
    }

    int single_dihedral_print = 5;                          //Number of single dihedral eigenvalues to print
    for(int i=0; i<N; i++){
        if(single_dihedral_basis_order[i] < single_dihedral_print){
            single_dihedral_print = single_dihedral_basis_order[i];
        }
    }

    cout << "=================================================================" << endl;
    cout << "          SINGLE ROTOR CALCULATION - CONFIGURATION DATA" << endl;
    cout << "=================================================================" << endl << endl;
    cout << "Number of rotors: " << N_rot << " (" << N << " dihedrals)" << endl;
    cout << "Max number of integration points: " << npt_int_single << endl;
    cout << "GSL QAG key: " << key_int_single << endl;
    cout << "Max integral errors:" << endl;
    cout << " -> Absolute: " << abs_int_single << endl;
    cout << " -> Relative: " << rel_int_single << endl << endl;
    cout << "DIHEDRAL ANGLES TABLE:" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << "| " << setw(3) << "i)" << " |" << setw(6) << "mins" << " |" <<
        setw(6) << "N Four" << " |" << setw(12) << "Barrier" << " |" <<
        setw(12) << "D(i-1)" << " |" << setw(12) << "D(i)" << " |" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    for(int i=0; i<N; i++){
        string label = to_string(i+1) + ")";
        cout << "| " << setw(3) << label << " |" << setw(6) << dihedral_num_mins[i] <<
        " |" << setw(6) << single_dihedral_basis_order[i] << " |" << scientific << setprecision(4)
        << setw(12) << dihedral_barrier[i] << " |" << setw(12) << D[i] <<
         " |" << setw(12) << D[i+1] << " |" << endl;
    }
    cout << "-----------------------------------------------------------------" << endl << endl;
    
    rotors::ISOLATED_SOLVER * Isolated_Basis_Set = new rotors::ISOLATED_SOLVER [N];

    for(int i=0; i<N; i++){
        rotors::ISOLATED_SOLVER solver(dihedral_num_mins[i], dihedral_barrier[i], single_dihedral_basis_order[i], npt_int_single, abs_int_single, rel_int_single, key_int_single);
        solver.solve();
        cout << "DIHEDRAL " << i+1 << endl;
        cout << "-----------------------------------------------------------------" << endl;
        cout << scientific << setprecision(10);
        cout << setw(3) << "N" << setw(20) << "EVEN      " << setw(20) << "ODD       " << endl;
        cout << "-----------------------------------------------------------------" << endl;
        for(int j=0; j<single_dihedral_print; j++){
            cout << setw(3) << j << setw(20) << solver.get_subspace_eigenval(j, true) << setw(20) << solver.get_subspace_eigenval(j, false) << endl;
        }
        cout << "-----------------------------------------------------------------" << endl << endl;
        Isolated_Basis_Set[i] = solver;
        if(debug_key == true){
            string debug_filename = "iso_" + to_string(i) + ".txt";
            ofstream debug_file;
            debug_file.open(debug_filename);
            int N_even = (single_dihedral_basis_order[i]%2==0)? single_dihedral_basis_order[i]/2 : 1+((single_dihedral_basis_order[i]-1)/2);
            int N_odd = (single_dihedral_basis_order[i]%2==0)? single_dihedral_basis_order[i]/2 : (single_dihedral_basis_order[i]-1)/2;
            debug_file << "# EVEN SUBSPACE" << endl;
            for(int c=0; c<N_even; c++){
                debug_file << c;
                for(int r=0; r<N_even; r++){
                    debug_file << '\t' << scientific << setprecision(12) << solver.get_subspace_eigvect(r, c, true);
                }
                debug_file << endl;
            }
            debug_file << "# ODD SUBSPACE" << endl;
            for(int c=0; c<N_odd; c++){
                debug_file << c;
                for(int r=0; r<N_odd; r++){
                    debug_file << '\t' << scientific << setprecision(12) << solver.get_subspace_eigvect(r, c, false);
                }
                debug_file << endl;
            }
            debug_file.close();
        }
    }

    //COUPLED SOLVER SECTION
    cout << "=================================================================" << endl;
    cout << "          COUPLED ROTOR CALCULATION - CONFIGURATION DATA" << endl;
    cout << "=================================================================" << endl << endl;
    cout << "Number of rotors: " << N_rot << " (" << N << " dihedrals)" << endl;
    cout << "Max number of integration points: " << npt_int_system << endl;
    cout << "GSL QAG key: " << key_int_system << endl;
    cout << "Max integral errors:" << endl;
    cout << " -> Absolute: " << abs_int_system << endl;
    cout << " -> Relative: " << rel_int_system << endl;
    
    if(scan_flag==true || locked_scan_flag==true){

        ofstream scan_file;
        scan_file.open("scan_report.txt");

        if(start_scan>end_scan){
            int buffer = end_scan;
            end_scan = start_scan;
            start_scan = buffer;
            step_scan = abs(step_scan);
        }

        cout << "=================================================================" << endl;
        cout << "            COUPLED ROTOR CALCULATION - SCAN MODE" << endl;
        cout << "=================================================================" << endl << endl;
        string scan_mode = (scan_flag==true)? "single" : "locked";
        cout << "Scan mode: " << scan_mode;
        if(scan_flag==true){
            cout << ", Selected dihedral: " << scan_index;
        }
        cout << endl;
        cout << "Start basis: " << start_scan << ", End basis: " << end_scan << ", Scan: " << step_scan << endl;

        int * scan_cutoff = new int [N];
        int scan_counter = 0;
        for(int basis_num=start_scan; basis_num<=end_scan; basis_num=basis_num+step_scan){
            scan_counter++;

            if(scan_flag==true){
                math_utils::copy_array<int>(composite_basis_set_cutoff, scan_cutoff, N);
                scan_cutoff[scan_index] = basis_num;
            }
            else{
                for(int i=0; i<N; i++){
                    scan_cutoff[i] = basis_num;
                }
            }

            M = 1;
            for(int i=0; i < N; i++){
                if(scan_cutoff[i] >= int(INT_MAX/M)){
                    cout << "ERROR (main.cpp): Integer overflow in generating the total number of basis functions (SCAN)" << endl;
                    exit(EXIT_FAILURE);
                }
                M *= scan_cutoff[i];
            }
            int coupled_dihedral_print = (M<10)? int(M/2)-1 : 5;    //Number of coupled dihedrals eigenvalues to print

            rotors::COUPLED_SOLVER System_Solver(N, D, scan_cutoff, Isolated_Basis_Set, npt_int_system, abs_int_system, rel_int_system, key_int_system);

            if(vqe_key==true){
                System_Solver.solve(true, false);
            }
            else{
                System_Solver.solve();
            }    
            
            cout << "=================================================================" << endl;
            cout << "                          SCAN STEP " << scan_counter << endl;
            cout << "=================================================================" << endl;
            cout << "Scan basis set dimension: " << basis_num << endl;
            cout << "Number of composite basis functions: " << M << endl;
            cout << "DIHEDRAL ANGLES TABLE:" << endl;
            cout << "-----------------------------------------------------------------" << endl;
            cout << "| " << setw(3) << "i)" << " |" << setw(6) << "mins" << " |" <<
                setw(6) << "N set" << " |" << setw(12) << "Barrier" << " |" <<
                setw(12) << "D(i-1)" << " |" << setw(12) << "D(i)" << " |" << endl;
            cout << "-----------------------------------------------------------------" << endl;
            for(int i=0; i<N; i++){
                string label = to_string(i+1) + ")";
                cout << "| " << setw(3) << label << " |" << setw(6) << dihedral_num_mins[i] <<
                " |" << setw(6) << scan_cutoff[i] << " |" << scientific << setprecision(4)
                << setw(12) << dihedral_barrier[i] << " |" << setw(12) << D[i] <<
                " |" << setw(12) << D[i+1] << " |" << endl;
            }
            cout << "-----------------------------------------------------------------" << endl;
            cout << "                  *** SYSTEM EIGENVALUES *** "<< endl;
            cout << "-----------------------------------------------------------------" << endl;
            cout << scientific << setprecision(10);
            cout << setw(3) << "N" << setw(20) << "GERADE     " << setw(20) << "UNGERADE    " << endl;
            cout << "-----------------------------------------------------------------" << endl;
            for(int j=0; j<coupled_dihedral_print; j++){
                cout << setw(3) << j << setw(20) << System_Solver.get_subspace_eigenval(j, true) << setw(20) << System_Solver.get_subspace_eigenval(j, false) << endl;
            }
            cout << endl << endl;
            
            scan_file << scan_counter << '\t' << basis_num << '\t' << scientific << setprecision(14) << System_Solver.get_subspace_eigenval(0, false) << endl;

            if(vqe_key==true){
                string vqe_filename = base_dir + "/VQE_" + to_string(basis_num) + ".txt";
                System_Solver.export_vqe_integrals(vqe_filename);
            }

            if(eigval_save_key==true){
                string eigval_filesname = base_dir + "/eigval_list_" + to_string(basis_num) + ".txt";
                System_Solver.export_eigenval_list(eigval_filesname);
            }

        }
        scan_file.close();
        delete[] scan_cutoff;
    }
    else{

        cout << "Number of composite basis functions: " << M << endl << endl;

        cout << "DIHEDRAL ANGLES TABLE:" << endl;
        cout << "-----------------------------------------------------------------" << endl;
        cout << "| " << setw(3) << "i)" << " |" << setw(6) << "mins" << " |" <<
            setw(6) << "N set" << " |" << setw(12) << "Barrier" << " |" <<
            setw(12) << "D(i-1)" << " |" << setw(12) << "D(i)" << " |" << endl;
        cout << "-----------------------------------------------------------------" << endl;
        for(int i=0; i<N; i++){
            string label = to_string(i+1) + ")";
            cout << "| " << setw(3) << label << " |" << setw(6) << dihedral_num_mins[i] <<
            " |" << setw(6) << composite_basis_set_cutoff[i] << " |" << scientific << setprecision(4)
            << setw(12) << dihedral_barrier[i] << " |" << setw(12) << D[i] <<
            " |" << setw(12) << D[i+1] << " |" << endl;
        }
        cout << "-----------------------------------------------------------------" << endl << endl;

        int coupled_dihedral_print = (M<10)? int(M/2)-1 : 5;    //Number of coupled dihedrals eigenvalues to print

        rotors::COUPLED_SOLVER System_Solver(N, D, composite_basis_set_cutoff, Isolated_Basis_Set, npt_int_system, abs_int_system, rel_int_system, key_int_system);

        if(vqe_key==true){
            System_Solver.solve(true, false);
        }
        else{
            System_Solver.solve();
        }    
        
        cout << "SYSTEM EIGENVALUES:" << endl;
        cout << "-----------------------------------------------------------------" << endl;
        cout << scientific << setprecision(10);
        cout << setw(3) << "N" << setw(20) << "GERADE     " << setw(20) << "UNGERADE    " << endl;
        cout << "-----------------------------------------------------------------" << endl;
        for(int j=0; j<coupled_dihedral_print; j++){
            cout << setw(3) << j << setw(20) << System_Solver.get_subspace_eigenval(j, true) << setw(20) << System_Solver.get_subspace_eigenval(j, false) << endl;
        }
        cout << "-----------------------------------------------------------------" << endl << endl;
        
        if(vqe_key==true){
            string vqe_filename = base_dir + "/VQE.txt";
            System_Solver.export_vqe_integrals(vqe_filename);
        }
        
        if(eigval_save_key==true){
            string eigval_filesname = base_dir + "/eigval_list.txt";
            System_Solver.export_eigenval_list(eigval_filesname);
        }

        if(debug_key == true){
            ofstream debug_file("comp.txt");
            int N_gerade = System_Solver.get_number_of_comp_functions(true);
            int N_ungerade = System_Solver.get_number_of_comp_functions(false);
            debug_file << "# GERADE COMPOSITE BASIS FUNCTIONS" << endl;
            for(int i=0; i<N_gerade; i++){
                debug_file << i;
                int *order_list = new int [N];
                bool *parity_list = new bool [N];
                System_Solver.get_base_function_composition(i, true, order_list, parity_list);
                for(int d=0; d<N; d++){
                    debug_file << '\t' << order_list[d];
                }
                for(int d=0; d<N; d++){
                    int parity_index = (parity_list[d] == true)? 1 : -1;
                    debug_file << '\t' << parity_index;
                }
                debug_file << '\n';
                delete[] order_list;
                delete[] parity_list;
            }
            debug_file << "# UNGERADE COMPOSITE BASIS FUNCTIONS" << endl;
            for(int i=0; i<N_ungerade; i++){
                debug_file << i;
                int *order_list = new int [N];
                bool *parity_list = new bool [N];
                System_Solver.get_base_function_composition(i, false, order_list, parity_list);
                for(int d=0; d<N; d++){
                    debug_file << '\t' << order_list[d];
                }
                for(int d=0; d<N; d++){
                    int parity_index = (parity_list[d] == true)? 1 : -1;
                    debug_file << '\t' << parity_index;
                }
                debug_file << '\n';
                delete[] order_list;
                delete[] parity_list;
            }
            debug_file.close();
        }
    }

    delete[] D;
    delete[] dihedral_barrier;
    delete[] dihedral_num_mins;
    delete[] single_dihedral_basis_order;
    delete[] composite_basis_set_cutoff;
    delete[] Isolated_Basis_Set;
    
    return 0;
}
