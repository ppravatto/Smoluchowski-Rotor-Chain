#ifndef ROTORS_H
#define ROTORS_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include "math_utils.h"
#include <omp.h>

namespace rotors{

    double fourier_basis(double theta, int n, bool parity){
        if(n==0) return 1./std::sqrt(2.*M_PI);
        double arg = (parity==true)? std::cos(n*theta) : std::sin(n*theta);
        return arg/std::sqrt(M_PI);
    }

    double fourier_basis_first_deriv(double theta, int n, bool parity){
        if(n==0) return 0.;
        double arg = (parity==true)? -n*std::sin(n*theta) : n*std::cos(n*theta);
        return arg/std::sqrt(M_PI);
    }

    double fourier_basis_second_deriv(double theta, int n, bool parity){
        if(n==0) return 0.;
        double arg = (parity==true)? -std::pow(n, 2.)*std::cos(n*theta) : -std::pow(n, 2.)*std::sin(n*theta);
        return arg/std::sqrt(M_PI);
    }

    class POTENTIAL{
        private:
            double barrier;
            int nmins;
        public:
            POTENTIAL(double barrier_, int nmins_){
                if(nmins_<=0 || nmins_>2){
                    std::cout << "ERROR (isolated_rotors.h - POTENTIAL): Invalid number of minima (" << nmins_ << ")" << std::endl;
                    exit(EXIT_FAILURE);
                }
                barrier = barrier_; nmins = nmins_;
            }

            double function(double theta){
                double arg = (nmins==1)? 1.-std::cos(theta) : std::cos(2.*theta)+1.;
                return barrier*arg/2.;
            }

            double first_deriv(double theta){
                double arg = (nmins==1)? std::sin(theta) : -2*std::sin(2.*theta);
                return barrier*arg/2.;
            }

            double second_deriv(double theta){
                double arg = (nmins==1)? std::cos(theta) : -4*std::cos(2.*theta);
                return barrier*arg/2.;
            }
    };

    struct carrier{
        int bra, ket;
        bool parity;
        POTENTIAL * MyPotential;
    };

    double matrix_element_integrand(double theta, void * pvoid){
        carrier p = * (carrier *) pvoid;
        POTENTIAL U = *(p.MyPotential);
        double arg = -fourier_basis_second_deriv(theta, p.ket, p.parity);
        arg += -0.5 * fourier_basis(theta, p.ket, p.parity) * U.second_deriv(theta);
        arg += 0.25 * fourier_basis(theta, p.ket, p.parity) * std::pow(U.first_deriv(theta), 2.);
        return fourier_basis(theta, p.bra, p.parity) * arg;
    }

    //Pre-declare functions to be friends of the ISOLATED_SOLVER class
    double one_rotor_term_integrand(double theta, void* pvoid);
    double two_rotors_deriv_term(double theta, void* pvoid);
    double two_rotor_poten_term(double theta, void* pvoid);

    class ISOLATED_SOLVER{
        private:
            bool solve_flag, initialization_flag;                               //Flag to verify the status of the diagonalization procedure
            int N, N_even, N_odd, nmins, key, npt;                              //Total number of basis functions, number of even and odd functions, number of potential minima, key for the QAG integrator and number of integration points
            double barrier, abs, rel;                                           //Barrier height, absolute and relative error for the integrator
            double *even_eigvect, *even_eigval, *odd_eigvect, *odd_eigval;      //Pointers to allocated memory

        protected:
            void check_init_error(){
                if(initialization_flag == false){
                    std::cout << "ERROR (isolated_rotors.h - ISOLATED_SOLVER): Object constructed by default constructor but not yet initialized" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }

            void check_access_errors(int order, int N_limit){
                if(solve_flag == false){
                    std::cout << "ERROR (isolated_rotors.h - ISOLATED_SOLVER): Eigenvalue problem not yet solved" << std::endl;
                    exit(EXIT_FAILURE);
                }
                if(order<0 || order>=N_limit){
                    std::cout << "SEGMENTATION FAULT (isolated_rotors.h - ISOLATED_SOLVER): Eigenvalue/Eigenfunction index out of bounds (core dumped)" << std::endl;
                    abort();
                }
            }

            void check_access_errors(int order, bool parity){
                int N_limit = (parity==true)? N_even : N_odd;
                check_access_errors(order, N_limit);
            }

            void transfer_data(const ISOLATED_SOLVER& s){
                if(s.initialization_flag == false){
                    std::cout << "ERROR (isolated_rotors.h - ISOLATED_SOLVER): Source class not yet inizialized" << std::endl;
                    exit(EXIT_FAILURE);
                }
                nmins = s.nmins; barrier = s.barrier; N = s.N; npt = s.npt; abs = s.abs;
                rel = s.rel; key = s.key; N_even = s.N_even; N_odd = s.N_odd; solve_flag = s.solve_flag;
                if(initialization_flag==false){
                    try{
                        even_eigval = new double [N_even];
                        even_eigvect = new double [N_even*N_even];
                        odd_eigval = new double [N_odd];
                        odd_eigvect = new double [N_odd*N_odd];
                        
                    }
                    catch(std::bad_alloc&){
                        std::cout << "ERROR (isolated_rotors.h - ISOLATED_SOLVER): Failure in allocating memory buffer" << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
                math_utils::copy_array<double>(s.even_eigval, even_eigval, N_even);
                math_utils::copy_array<double>(s.even_eigvect, even_eigvect, N_even*N_even);
                math_utils::copy_array<double>(s.odd_eigval, odd_eigval, N_odd);
                math_utils::copy_array<double>(s.odd_eigvect, odd_eigvect, N_odd*N_odd);
            }

        public:
            ISOLATED_SOLVER(){
                initialization_flag = false;
            }

            ISOLATED_SOLVER(int nmins_, double barrier_, int N_, int npt_, double abs_, double rel_, int key_){
                nmins = nmins_; barrier = barrier_; N = N_; npt = npt_; abs = abs_; rel = rel_; key = key_;
                N_even = (N%2==0)? N/2 : 1+((N-1)/2);
                N_odd = (N%2==0)? N/2 : (N-1)/2;
                try{
                    even_eigval = new double [N_even];
                    even_eigvect = new double [N_even*N_even];
                    odd_eigval = new double [N_odd];
                    odd_eigvect = new double [N_odd*N_odd];
                    
                }
                catch(std::bad_alloc&){
                    std::cout << "ERROR (isolated_rotors.h - ISOLATED_SOLVER): Failure in allocating memory buffer" << std::endl;
                    exit(EXIT_FAILURE);
                }
                solve_flag = false;
                initialization_flag = true;
            }

            ISOLATED_SOLVER(const ISOLATED_SOLVER& s){
                initialization_flag = false;
                transfer_data(s);
                initialization_flag = s.initialization_flag;
            }
            
            ISOLATED_SOLVER& operator = (const ISOLATED_SOLVER& s){
                transfer_data(s);
                initialization_flag = true;
                return *this;
            }

            ~ISOLATED_SOLVER(){
                if(initialization_flag == true){
                    delete[] even_eigvect;
                    delete[] even_eigval;
                    delete[] odd_eigvect;
                    delete[] odd_eigval;
                }
            }

            double get_matrix_element(int row, int col, bool parity){
                check_init_error();
                POTENTIAL U(barrier, nmins);
                carrier data;
                data.bra = (parity==true)? row : row+1;
                data.ket = (parity==true)? col : col+1;;
                data.parity = parity; data.MyPotential = &U;
                return math_utils::QAG_integrator(&matrix_element_integrand, &data, -M_PI, M_PI, npt, abs, rel, key);
            }

            void solve(){
                check_init_error();
                #ifdef _OPENMP
                    #pragma omp parallel for collapse(2) schedule(dynamic)
                    for(int col=0; col<N_even; col++){
                        for(int row=0; row<N_even; row++){
                            even_eigvect[row+N_even*col] = get_matrix_element(row, col, true);
                        }
                    }
                    #pragma omp parallel for collapse(2) schedule(dynamic)
                    for(int col=0; col<N_odd; col++){
                        for(int row=0; row<N_odd; row++){
                            odd_eigvect[row+N_odd*col] = get_matrix_element(row, col, false);
                        }
                    }
                #else
                    for(int col=0; col<N_even; col++){
                        for(int row=0; row<N_even; row++){
                            even_eigvect[row+N_even*col] = get_matrix_element(row, col, true);
                        }
                    }
                    for(int col=0; col<N_odd; col++){
                        for(int row=0; row<N_odd; row++){
                            odd_eigvect[row+N_odd*col] = get_matrix_element(row, col, false);
                        }
                    }
                #endif
                math_utils::matrix_symmetrize(even_eigvect, N_even);
                math_utils::matrix_symmetrize(odd_eigvect, N_odd);
                math_utils::symm_real_eigsys(even_eigvect, N_even, even_eigval);
                math_utils::symm_real_eigsys(odd_eigvect, N_odd, odd_eigval);

                solve_flag = true;
            }

            double get_subspace_eigvect(int element, int order, bool parity){
                check_access_errors(order, parity);
                return (parity==true)? even_eigvect[element+N_even*order] : odd_eigvect[element+N_odd*order];
            }

            double get_subspace_eigfunc_deriv(double theta, int order, bool parity, int deriv_order){
                check_init_error();
                check_access_errors(order, parity);
                double value = 0.;
                int N_limit = (parity==true)? N_even : N_odd;
                for(int row=0; row<N_limit; row++){
                    int index = (parity==true)? row : row+1;
                    double eigvect = (parity==true)? even_eigvect[row+N_even*order] : odd_eigvect[row+N_odd*order];
                    if(deriv_order==0){
                        value += eigvect*fourier_basis(theta, index, parity);
                    }
                    else if(deriv_order==1){
                        value += eigvect*fourier_basis_first_deriv(theta, index, parity);
                    }
                    else if(deriv_order==2){
                        value += eigvect*fourier_basis_second_deriv(theta, index, parity);
                    }
                    else{
                        std::cout << "ERROR (get_subspace_eigfunc_deriv.h - ISOLATED_SOLVER): Derivative order not available" << std::endl;
                    }
                }
                return value;
            }

            double get_subspace_eigfunc(double theta, int order, bool parity){
                return get_subspace_eigfunc_deriv(theta, order, parity, 0);
            }

            double get_subspace_eigenval(int order, bool parity){
                check_init_error();
                check_access_errors(order, parity);
                return (parity==true)? even_eigval[order] : odd_eigval[order];
            }
            
            friend double one_rotor_term_integrand(double, void*);
            friend double two_rotors_deriv_term(double, void*);
            friend double two_rotor_poten_term(double, void*);
    };

    struct coupled_carrier{
        int index, bra, ket, N;
        bool bra_parity, ket_parity;
        double * D;
        ISOLATED_SOLVER * IsolatedBasis;
    };

    double one_rotor_term_integrand(double theta, void* pvoid){
        coupled_carrier p = *(coupled_carrier *) pvoid;
        ISOLATED_SOLVER Dihedral = p.IsolatedBasis[p.index];
        POTENTIAL U(Dihedral.barrier, Dihedral.nmins);
        double arg = (0.5*U.second_deriv(theta) - 0.25*std::pow(U.first_deriv(theta), 2.)) * Dihedral.get_subspace_eigfunc(theta, p.ket, p.ket_parity);
        arg += Dihedral.get_subspace_eigfunc_deriv(theta, p.ket, p.ket_parity, 2);
        double coeff = p.D[p.index]+p.D[p.index+1];
        return -coeff*Dihedral.get_subspace_eigfunc(theta, p.bra, p.bra_parity)*arg;
    }
    
    double two_rotors_deriv_term(double theta, void* pvoid){
        coupled_carrier p = *(coupled_carrier *) pvoid;
        ISOLATED_SOLVER Dihedral = p.IsolatedBasis[p.index];
        return Dihedral.get_subspace_eigfunc(theta, p.bra, p.bra_parity)*Dihedral.get_subspace_eigfunc_deriv(theta, p.ket, p.ket_parity, 1);
    }

    double two_rotor_poten_term(double theta, void* pvoid){
        coupled_carrier p = *(coupled_carrier *) pvoid;
        ISOLATED_SOLVER Dihedral = p.IsolatedBasis[p.index];
        POTENTIAL U(Dihedral.barrier, Dihedral.nmins);
        return Dihedral.get_subspace_eigfunc(theta, p.bra, p.bra_parity)*U.first_deriv(theta)*Dihedral.get_subspace_eigfunc(theta, p.ket, p.ket_parity);
    }

    class COUPLED_SOLVER{
        private:
            int N, N_comp, N_gerade, N_ungerade;                    //Number of dihedral angles, total number of composite functions, number of gerade and ungerade composite functions;
            int npt, key;                                           //Number of integration points and key for the QAG integrator
            double abs, rel;                                        //Absolute and relative errors for the integrator
            bool symmetry_flag, parity, solve_flag;                 //Flag to activate the symmetry constrained active space, flag to select the parity of the symmetry constrained active space, flag to verify the status of the diagonalization procedure
            bool vqe_status, vqe_parity;                            //Flag to indicate the initialization status of the Hamiltonian matrix, Flag to enable VQE copy of the data, VQE parity selection
            int *N_cutoff, *gerade_list, *ungerade_list;            //Pointer to the number of single dihedral eigenfunction to be included for angle (maximum order), list of composite functions order of a given inversion symmetry
            double *D;                                              //Pointer to the diffusion coefficient list
            double *eigval_gerade, *eigvect_gerade;                 //Pointer to the eigenvalues and eigenvectors arrays for the grade space
            double *eigval_ungerade, *eigvect_ungerade;             //Pointer to the eigenvalues and eigenvectors arrays for the ungrade space
            double *vqe_matrix;
            ISOLATED_SOLVER *IsolatedBasis;                         //Pointer to the object defining the basis set for the isolated rotor

            void check_access_errors(int order, int N_limit){
                if(order<0 || order>=N_limit){
                    std::cout << "SEGMENTATION FAULT (coupled_rotors.h - COUPLED_SOLVER): Basis function index out of bounds (core dumped)" << std::endl;
                    abort();
                }
            }

            void check_solve_status(){
                if(solve_flag==false){
                    std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER):  Eigenvalue problem not yet solved" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }

        protected:
            void order_decomposition(int order, int* order_list, bool* parity_list){
                int value = order;
                for(int i=N-1; i>=0; i--){
                    int div = 1;
                    for(int j=0; j<i; j++){
                        div *= N_cutoff[j];
                    }
                    int rem = value%div;
                    int index = (value-rem)/div;
                    parity_list[i] = (index%2==0)? true : false; 
                    order_list[i] = (index-(index%2))/2;
                    value = rem;
                }
            }

            bool get_parity(int order){
                int* temp_order_list; 
                bool* temp_parity_list;
                try{
                    temp_order_list = new int [N];
                    temp_parity_list = new bool [N];
                }
                catch(std::bad_alloc&){
                    std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER): Failure in allocating temporary memory buffer" << std::endl;
                    exit(EXIT_FAILURE);
                }
                order_decomposition(order, temp_order_list, temp_parity_list);
                int func_parity = 1;
                for(int j=0; j<N; j++){
                    func_parity *= (temp_parity_list[j]==true)? 1 : -1;
                }
                delete[] temp_order_list;
                delete[] temp_parity_list;
                return (func_parity>0)? true : false;
            }

            void constructor(int N_, double* D_, int* N_cutoff_, ISOLATED_SOLVER* IsolatedBasis_, int npt_, double abs_, double rel_, int key_){
                N = N_; D = D_; N_cutoff = N_cutoff_; IsolatedBasis = IsolatedBasis_; npt = npt_; abs = abs_; rel = rel_; key = key_;
                vqe_status = false;
                N_comp = 1; N_gerade = 0; N_ungerade=0;
                for(int i=0; i<N; i++){
                    N_comp *= N_cutoff[i];
                }
                for(int i=0; i<N_comp; i++){
                    bool func_parity = get_parity(i);
                    if(func_parity==true){
                        N_gerade++;
                    }
                    else{
                        N_ungerade++;
                    }
                }
                try{
                    gerade_list = new int [N_gerade];
                    ungerade_list = new int [N_ungerade];
                }
                catch(std::bad_alloc&){
                    std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER): Failure in allocating memory buffer" << std::endl;
                    exit(EXIT_FAILURE);
                }
                int index_gerade=0, index_ungerade=0;
                for(int i=0; i<N_comp; i++){
                    bool func_parity = get_parity(i);
                    if(func_parity == true){
                        gerade_list[index_gerade] = i;
                        index_gerade++;
                    }
                    else{
                        ungerade_list[index_ungerade] = i;
                        index_ungerade++;
                    }
                }
                try{
                    if(symmetry_flag==false || (symmetry_flag==true && parity==true)){
                        eigval_gerade = new double [N_gerade];
                        eigvect_gerade = new double [N_gerade*N_gerade];
                    }
                    if(symmetry_flag==false || (symmetry_flag==true && parity==false)){
                        eigval_ungerade = new double [N_ungerade];
                        eigvect_ungerade = new double [N_ungerade*N_ungerade];
                    }
                }
                catch(std::bad_alloc&){
                    std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER): Failure in allocating memory buffer" << std::endl;
                    exit(EXIT_FAILURE);
                }
            }

        public:
            COUPLED_SOLVER(int N_, double* D_, int* N_cutoff_, ISOLATED_SOLVER* IsolatedBasis_, int npt_, double abs_, double rel_, int key_){
                symmetry_flag = false;
                constructor(N_, D_, N_cutoff_, IsolatedBasis_, npt_, abs_, rel_, key_);
            }

            COUPLED_SOLVER(int N_, double* D_, int* N_cutoff_, ISOLATED_SOLVER* IsolatedBasis_, int npt_, double abs_, double rel_, int key_, bool parity_){
                symmetry_flag = true;
                parity = parity_;
                constructor(N_, D_, N_cutoff_, IsolatedBasis_, npt_, abs_, rel_, key_);
            }

            ~COUPLED_SOLVER(){
                delete[] gerade_list;
                delete[] ungerade_list;
                if(symmetry_flag==false || (symmetry_flag==true && parity==true)){
                    delete[] eigval_gerade;
                    delete[] eigvect_gerade;
                }
                if(symmetry_flag==false || (symmetry_flag==true && parity==false)){
                    delete[] eigval_ungerade;
                    delete[] eigvect_ungerade;
                }
                if(vqe_status==true) delete[] vqe_matrix;
            }

            double comp_basis_func(double* theta, int order){
                check_access_errors(order, N_comp);
                int* one_rot_order;
                bool* one_rot_parity;
                try{
                    one_rot_order = new int [N];
                    one_rot_parity = new bool [N];
                }
                catch(std::bad_alloc&){
                    std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER): Failure in allocating temporary memory buffer" << std::endl;
                    exit(EXIT_FAILURE);
                }
                order_decomposition(order, one_rot_order, one_rot_parity);
                double func = 1.;
                for(int i=0; i<N; i++){
                    func *= (IsolatedBasis[i]).get_subspace_eigfunc(theta[i], one_rot_order[i], one_rot_parity[i]);
                }
                delete[] one_rot_order; 
                delete[] one_rot_parity;
                return func;
            }

            double comp_basis_func(double* theta, int order, bool function_parity){
                int N_limit = (function_parity==true)? N_gerade : N_ungerade;
                check_access_errors(order, N_limit);
                int effective_order = (function_parity==true)? gerade_list[order] : ungerade_list[order];
                return comp_basis_func(theta, effective_order);
            }

            bool comp_basis_func_parity(int order){
                check_access_errors(order, N_comp);
                return get_parity(order);
            }

            int get_number_of_comp_functions(bool parity){
                return (parity==true)? N_gerade : N_ungerade;
            }

            void get_base_function_composition(int index, bool parity, int* order_list, bool* parity_list){
                int N_limit = (parity==true)? N_gerade : N_ungerade;
                check_access_errors(index, N_limit);
                int effective_index = (parity==true)? gerade_list[index] : ungerade_list[index];
                order_decomposition(effective_index, order_list, parity_list);
            }

            double get_matrix_element(int row, int col, bool element_parity){
                int N_limit = (element_parity==true)? N_gerade : N_ungerade;
                check_access_errors(row, N_limit);
                check_access_errors(col, N_limit);
                int effective_row = (element_parity==true)? gerade_list[row] : ungerade_list[row];
                int effective_col = (element_parity==true)? gerade_list[col] : ungerade_list[col];
                int *order_list_row, *order_list_col;
                bool *parity_list_row, *parity_list_col;
                try{
                    order_list_row = new int [N];
                    order_list_col = new int [N];
                    parity_list_row = new bool [N];
                    parity_list_col = new bool [N];
                }
                catch(std::bad_alloc&){
                    std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER): Failure in allocating temporary memory buffer" << std::endl;
                    exit(EXIT_FAILURE);
                }
                order_decomposition(effective_row, order_list_row, parity_list_row);
                order_decomposition(effective_col, order_list_col, parity_list_col);
                coupled_carrier data;
                data.D = D; data.N = N; data.IsolatedBasis = IsolatedBasis;
                double sum = 0.;
                for(int n=0; n<N; n++){
                    data.index = n;
                    data.bra = order_list_row[n]; data.bra_parity = parity_list_row[n];
                    data.ket = order_list_col[n]; data.ket_parity = parity_list_col[n];
                    
                    bool ortogonality_flag_single = false;
                    bool ortogonality_flag_double = false;
                    for(int i=0; i<N; i++){
                        if((order_list_row[i] != order_list_col[i]) || (parity_list_row[i] != parity_list_col[i])){
                            ortogonality_flag_single = true;
                            if(i!=n && i!=n+1) ortogonality_flag_double = true;
                        }
                    }

                    if(ortogonality_flag_single==false){
                        sum += math_utils::QAG_integrator(&one_rotor_term_integrand, &data, -M_PI, M_PI, npt, abs, rel, key);
                    }
                    
                    if(n<N-1 && ortogonality_flag_double==false){
                        double deriv_term_1 = math_utils::QAG_integrator(&two_rotors_deriv_term, &data, -M_PI, M_PI, npt, abs, rel, key);
                        double pot_term_1 = math_utils::QAG_integrator(&two_rotor_poten_term, &data, -M_PI, M_PI, npt, abs, rel, key);
                        data.index = n+1;
                        data.bra = order_list_row[n+1]; data.bra_parity = parity_list_row[n+1];
                        data.ket = order_list_col[n+1]; data.ket_parity = parity_list_col[n+1];
                        double deriv_term_2 = math_utils::QAG_integrator(&two_rotors_deriv_term, &data, -M_PI, M_PI, npt, abs, rel, key);
                        double pot_term_2 = math_utils::QAG_integrator(&two_rotor_poten_term, &data, -M_PI, M_PI, npt, abs, rel, key);
                        sum += 2.*D[n+1]*(deriv_term_1*deriv_term_2 - 0.25*pot_term_1*pot_term_2);
                        
                    }
                    
                }

                delete[] order_list_row;
                delete[] parity_list_row;
                delete[] order_list_col;
                delete[] parity_list_col;
                return sum;
            }

            void solve(bool vqe_, bool vqe_parity_){
                if(vqe_==true && symmetry_flag==true && parity != vqe_parity_){
                    std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER): The selected VQE parity does not match the class constructor" << std::endl;
                    exit(EXIT_FAILURE);
                }

                if(vqe_==true){
                    vqe_parity = vqe_parity_;
                    try{
                        if(vqe_parity ==false){
                            vqe_matrix = new double [N_gerade*N_gerade];
                        }
                        else{
                            vqe_matrix = new double [N_ungerade*N_ungerade];
                        }
                    }
                    catch(std::bad_alloc&){
                        std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER): Failure in allocating temporary memory buffer for VQE" << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }

                if(symmetry_flag==false || (symmetry_flag==true && parity==true)){
                    #ifdef _OPENMP
                        #pragma omp parallel for collapse(2) schedule(dynamic)
                        for(int row=0; row<N_gerade; row++){
                            for(int col=0; col<N_gerade; col++){
                                eigvect_gerade[row+N_gerade*col] = get_matrix_element(row, col, true);
                            }
                        }
                    #else
                        for(int row=0; row<N_gerade; row++){
                            for(int col=0; col<N_gerade; col++){
                                eigvect_gerade[row+N_gerade*col] = get_matrix_element(row, col, true);
                            }
                        }
                    #endif
                    math_utils::matrix_symmetrize(eigvect_gerade, N_gerade);
                    if(vqe_==true && vqe_parity==true) math_utils::copy_array<double>(eigvect_gerade, vqe_matrix, N_gerade*N_gerade);
                    math_utils::symm_real_eigsys(eigvect_gerade, N_gerade, eigval_gerade);
                }

                if(symmetry_flag==false || (symmetry_flag==true && parity==false)){
                    #ifdef _OPENMP
                        #pragma omp parallel for collapse(2) schedule(dynamic)
                        for(int row=0; row<N_ungerade; row++){
                            for(int col=0; col<N_ungerade; col++){
                                eigvect_ungerade[row+N_ungerade*col] = get_matrix_element(row, col, false);
                            }
                        }
                    #else
                        for(int row=0; row<N_ungerade; row++){
                            for(int col=0; col<N_ungerade; col++){
                                eigvect_ungerade[row+N_ungerade*col] = get_matrix_element(row, col, false);
                            }
                        }
                    #endif
                    math_utils::matrix_symmetrize(eigvect_ungerade, N_ungerade);
                    if(vqe_==true && vqe_parity==false) math_utils::copy_array<double>(eigvect_ungerade, vqe_matrix, N_ungerade*N_ungerade);
                    math_utils::symm_real_eigsys(eigvect_ungerade, N_ungerade, eigval_ungerade);
                }

                vqe_status = vqe_;
                solve_flag = true;
            }

            void solve(){
                solve(false, false);
            }

            double get_subspace_eigfunc(double* theta, int order, bool function_parity){
                check_solve_status();
                if(symmetry_flag==true && (function_parity != parity)){
                    std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER): The required eigenfunction has not the selected symmetry" << std::endl;
                    exit(EXIT_FAILURE);
                }
                int N_limit = (function_parity==true)? N_gerade : N_ungerade;
                check_access_errors(order, N_limit);
                double sum = 0.;
                for(int row=0; row<N_limit; row++){
                    double eigvect = (function_parity==true)? eigvect_gerade[row+N_gerade*order] : eigvect_ungerade[row+N_ungerade*order];
                    sum += eigvect*comp_basis_func(theta, row, function_parity);
                }
                return sum;
            }

            double get_subspace_eigenval(int order, bool function_parity){
                check_solve_status();
                if(symmetry_flag==true && (function_parity != parity)){
                    std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER): The required eigenvalue has not the selected symmetry" << std::endl;
                    exit(EXIT_FAILURE);
                }
                int N_limit = (function_parity==true)? N_gerade : N_ungerade;
                check_access_errors(order, N_limit);
                return (function_parity==true)? eigval_gerade[order] : eigval_ungerade[order];
            }

            void export_eigenval_list(std::string filename){
                std::ofstream datafile;
                datafile.open(filename);
                if(symmetry_flag==true){
                    int N_limit = (parity==true)? N_gerade : N_ungerade;
                    for(int i=0; i<N_limit; i++){
                        std::string parity_label = (parity==true)? "G" : "U";
                        double eigval = (parity==true)? eigval_gerade[i] : eigval_ungerade[i];
                        datafile << i << '\t' << parity_label << '\t' << eigval << std::endl;
                    }
                }
                else{
                    int index_gerade=0, index_ungerade=0, index_sum=0;
                    while(index_sum < N_comp){
                        if(index_gerade>=N_gerade || index_ungerade>=N_ungerade){
                            break;
                        }
                        if(eigval_gerade[index_gerade] < eigval_ungerade[index_ungerade]){
                            datafile << index_gerade << '\t' << "G" << '\t' << eigval_gerade[index_gerade] << std::endl;
                            index_gerade++;
                        }
                        else{
                            datafile << index_ungerade << '\t' << "U" << '\t' << eigval_ungerade[index_ungerade] << std::endl;
                            index_ungerade++;
                        }
                        index_sum = index_gerade + index_ungerade;
                    }
                    while(index_sum < N_comp){
                        if(index_gerade>=N_gerade){
                            datafile << index_ungerade << '\t' << "U" << '\t' << eigval_ungerade[index_ungerade] << std::endl;
                            index_ungerade++;
                        }
                        else{
                            datafile << index_gerade << '\t' << "G" << '\t' << eigval_gerade[index_gerade] << std::endl;
                            index_gerade++;
                        }
                        index_sum = index_gerade + index_ungerade;
                    }
                }
                datafile.close();
            }

            void export_vqe_integrals(std::string filename){
                if(vqe_status == false){
                    std::cout << "ERROR (coupled_rotors.h - COUPLED_SOLVER): VQE data not available" << std::endl;
                    exit(EXIT_FAILURE);
                }
                int N_limit = (vqe_parity==true)? N_gerade : N_ungerade;
                std::ofstream output;
                output.open(filename);
                output << std::scientific << std::setprecision(15);
                for(int row=0; row<N_limit; row++){
                    for(int col=0; col<N_limit; col++){
                        output << vqe_matrix[row+N_limit*col] << '\t';
                    }
                    output << std::endl;
                }
                output.close();
            }
    };
}

#endif