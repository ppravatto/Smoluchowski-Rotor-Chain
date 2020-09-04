#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include <string>
#include <limits.h>
#include <omp.h>
#include "rotors.h"
#include "input_parser.h"

using namespace std;

int main(int argc, char** argv){

    if(argc==1){
        std::cout << "ERROR (main.cpp): No input file selected" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::string filename = argv[1];
    std::cout << "Selected input file: " << filename << std::endl << std::endl;

    input::INPUT_PARSER datafile(filename);             //Load input file with the input parser
    int N_rot = datafile.load();                        //Set total number of rotors
    int N = N_rot-1;                                    //Number od dihedral angles

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
    int coupled_dihedral_print = (M<10)? int(M/2)-1 : 5;    //Number of coupled dihedrals eigenvalues to print
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
    cout << " -> Relative: " << rel_int_single << endl;
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
    
    rotors::COUPLED_SOLVER System_Solver(N, D, composite_basis_set_cutoff, Isolated_Basis_Set, npt_int_system, abs_int_system, rel_int_system, key_int_system);
    System_Solver.solve(true, false);
    
    cout << "SYSTEM EIGENVALUES:" << endl;
    cout << "-----------------------------------------------------------------" << endl;
    cout << scientific << setprecision(10);
    cout << setw(3) << "N" << setw(20) << "GERADE     " << setw(20) << "UNGERADE    " << endl;
    cout << "-----------------------------------------------------------------" << endl;
    for(int j=0; j<coupled_dihedral_print; j++){
        cout << setw(3) << j << setw(20) << System_Solver.get_subspace_eigenval(j, true) << setw(20) << System_Solver.get_subspace_eigenval(j, false) << endl;
    }
    cout << "-----------------------------------------------------------------" << endl << endl;
    
    System_Solver.export_vqe_integrals("VQE.txt");
    
    delete[] D;
    delete[] dihedral_barrier;
    delete[] dihedral_num_mins;
    delete[] single_dihedral_basis_order;
    delete[] composite_basis_set_cutoff;
    delete[] Isolated_Basis_Set;
    
    return 0;
}
