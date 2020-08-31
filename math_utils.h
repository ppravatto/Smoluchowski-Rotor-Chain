#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <iostream>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_integration.h>

extern "C"{
    extern int dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
}

namespace math_utils{

    //Routine to print a generic rectangular matrix of doubles
    void print_matrix(double *M, int n_rows, int n_cols){
        for(int row=0; row<n_rows; row++){
            for(int col=0; col<n_cols; col++){
                std::cout << std::scientific << std::setprecision(4) << std::setw(12) << M[row+n_rows*col];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    //Overload of the previous routine to print square matrices
    void print_matrix(double *M, int n){
        print_matrix(M, n, n);
    }

    //The routine diagonalizes a NxN square symmetric real matrix. The matrix is overwritten by the eigenvectors and the eigenvalues ar stored in "eigenval"
    void symm_real_eigsys(double *eigvect, int N, double* eigval){
        char jobz = 'V';
        char uplo = 'U';
        int lwork = 3*N;
        int info;
        double * work;
        try{
            work = new double [lwork];
        }
        catch(std::bad_alloc&){
            std::cout << "ERROR (math_utils.h - symm_real_eigsys): Failure in allocating LAPACK workspace memory buffer" << std::endl;
            exit(EXIT_FAILURE);
        }
        dsyev_(&jobz, &uplo, &N, eigvect, &N, eigval, work, &lwork, &info);
        delete[] work;
        if(info!=0){
            if(info>0){
                std::cout << "ERROR (math_utils.h - symm_real_eigsys): [DSYEV] algorithm did not converge" << std::endl;
            }
            else{
                std::cout << "ERROR (math_utils.h - symm_real_eigsys): [DSYEV] illegal value at argument " << -info << std::endl;
            }
            exit(EXIT_FAILURE);
        }
    }

    //The routine symmetrize a square matrix
    void matrix_symmetrize(double* M, int N){
        for(int col=0; col<N; col++){
            for(int row=0; row<col; row++){
                double average = 0.5*M[row+N*col] + 0.5*M[col+N*row];
                M[row+N*col]= average; M[col+N*row]= average;
            }
        }
    }

    //Template routine to copy an array given the starting pointer
    template <class MyType>
    void copy_array(MyType* source, MyType* target, int length){
        if(length<=0){
            std::cout << "ERROR (math_utils.h - copy_array): Illegal length argument" << std::endl;
            exit(EXIT_FAILURE);
        }
        for(int i=0; i<length; i++){
            target[i] = source[i];
        }
    }

    //The routine computes the integral of a double f(double, void*) in a confined interval x_min, x_max using the QAG algorithm
    double QAG_integrator(double (*f)(double, void*), void * pvoid, double x_min, double x_max, int npt, double abs, double rel, int key){
        double result, error;
        gsl_function integrand;
        integrand.function = f;
        integrand.params = pvoid;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(npt);
        gsl_integration_qag(&integrand, x_min, x_max, abs, rel, npt, key, w, &result, &error);
        gsl_integration_workspace_free(w);
        return result;
    }
    
}
#endif