#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <complex>
#include <cmath>
#include <cassert>
//#include "mkl_spblas.h"
//#include <mkl_types.h>
//#include <mkl_cblas.h>
//#include <mkl_lapacke.h>
//#include <Eigen/Dense>
//#include <Eigen/Eigenvalues>
#include "tensors.h"
#include "Parameters_BHZ.h"
#include "Hamiltonian_BHZ.h"
#include "Observables_BHZ.h"

using namespace std;
//using namespace Eigen;

int main(int argc, char *argv[]){

    string inputfile = argv[1];

    //Calling out the classes and their constructors:
    Parameters_BHZ Parameters_BHZ_;
    Parameters_BHZ_.Initialize(inputfile);
        
    Hamiltonian_BHZ Hamiltonian_BHZ_(Parameters_BHZ_);
    Hamiltonian_BHZ_.connectionMatrix();
    Hamiltonian_BHZ_.Diagonalizer();

    Observables_BHZ Observables_BHZ_(Parameters_BHZ_,Hamiltonian_BHZ_);
    if(Parameters_BHZ_.get_ldoe){
        Observables_BHZ_.calculateLDOE();
    }
    if(Parameters_BHZ_.get_dos){
        Observables_BHZ_.calculateDOS();
    }
    if(Parameters_BHZ_.get_sc){
        Observables_BHZ_.calculateSpinCurrents();
        Observables_BHZ_.calculateEdgeConductance();
    }
    if(Parameters_BHZ_.get_Akxw){
    //    Observables_BHZ_.calculateBmat();
        Observables_BHZ_.calculateAkxw();
    }
    if(Parameters_BHZ_.get_Akyw){
    //    Observables_BHZ_.calculateBmat();
        Observables_BHZ_.calculateAkyw();
    //    Observables_BHZ_.calculateAkywNew();
    }
    if(Parameters_BHZ_.get_ky_resolved_Akxw){
    //    Observables_BHZ_.calculateBmat();
    //    Observables_BHZ_.calculateYMomentumResolvedAkxw();
    }
    if(Parameters_BHZ_.get_wave_fn){
        Observables_BHZ_.calculateWaveFunctions();
    }
    if(Parameters_BHZ_.get_ldos){
        Observables_BHZ_.calculateLDOS();
    }

//    Observables_BHZ_.B_mat.clear();

    return 0;
}
