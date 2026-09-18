#include "tensors.h"
#include "Parameters_BHZ.h"
extern "C" {
//    #include <lapacke.h>
//    #include <cblas.h>
    void zheev_(char *, char *, int *, std::complex<double> *, int *, double *, std::complex<double> *, int *, double *, int *);
}
//using namespace Eigen;

#ifndef Hamiltonian_BHZ_class
#define Hamiltonian_BHZ_class
#define PI acos(-1.0)

class Hamiltonian_BHZ{

    public:
    //Constructor:
    Hamiltonian_BHZ(Parameters_BHZ &Parameters_BHZ__): Parameters_BHZ_(Parameters_BHZ__){
        Initialize();
    }
    void Initialize();

    void addOnsiteTerms(int r, int orb, int spin);
    void addHoppingXTerm(int r, int rpx, int orb);
    void addHoppingYTerm(int r, int rpy, int orb);
    void addOrbMixingXTerm(int r, int rpx, int orb, int spin);
    void addOrbMixingYTerm(int r, int rpy, int orb, int spin);

    void addTBC_HoppingYTerm(int r, int rpy, int orb, int spin);
    void addTBC_OrbMixingYTerm(int r, int rpy, int orb, int spin);
    void addPinningTerm(int r, int orb, int spin);
    
    void connectionMatrix();
    void Diagonalizer();
    int bar_function(int val);
    int makeIndex(int rx, int ry, int orb, int spin);

    void updownOccupiedSubspace();

    Parameters_BHZ &Parameters_BHZ_;
    
    int lx_, ly_, wx_;
    int N_cells_,N_orbs_,N_spin_;
    int N_particles_, H_size, half_size_;
    double A_, B_, M_, D_, C_, Vo_;

    Mat_1_doub evals_;
    Mat_2_Complex_doub C_mat;
    //Mat_1_Complex_doub evecs_;
    Mat_2_Complex_doub evecs_;
};

void Hamiltonian_BHZ::Initialize(){

    A_ = Parameters_BHZ_.A_val;
    B_ = Parameters_BHZ_.B_val;
    M_ = Parameters_BHZ_.M_val;
    D_ = Parameters_BHZ_.D_val;
    C_ = Parameters_BHZ_.C_val;

    lx_ = Parameters_BHZ_.Lx_;
    ly_ = Parameters_BHZ_.Ly_;
    wx_ = Parameters_BHZ_.W_;
    Vo_ = Parameters_BHZ_.Vo_val;
    
    N_cells_ = Parameters_BHZ_.Total_Cells;
    N_orbs_ = Parameters_BHZ_.N_Orbs;
    N_spin_ = 2;
    H_size = Parameters_BHZ_.Ham_Size;
    half_size_=(int) (H_size/2);

    evals_.resize(H_size);
    //evecs_.resize(H_size*H_size);
    evecs_.resize(H_size);
    C_mat.resize(H_size);
    for(int i=0;i<H_size;i++){
        evecs_[i].resize(H_size);
        C_mat[i].resize(H_size);
    }
/*
cout << "A: " << A_ << ", B: " << B_ << ", M: " << M_ << endl;
cout << "lx: " << lx_ << ", ly: " << ly_ << ", wx: " << wx_ << endl;
cout << "N_cells: " << N_cells_ << ", N_orbs: " << N_orbs_ << ", N_spin: " << N_spin_ << endl;
cout << "H_size: " << H_size << ", half_size: " << half_size_ << endl;
*/   
}

int Hamiltonian_BHZ::makeIndex(int rx, int ry, int orb, int spin){
	int r_val;
	r_val = spin * half_size_ + (orb + 2 * ry + 2 * ly_ * rx);
	return r_val;
}

int Hamiltonian_BHZ::bar_function(int val){
    int return_val;
    if (val == 0){      return_val = 1; }
	else if(val ==1){   return_val = 0; }

	return return_val;
}

void Hamiltonian_BHZ::addOnsiteTerms(int r, int orb, int spin){
	
//cout << "Adding onsite term: r=" << r << ", orb=" << orb << ", spin=" << spin << endl;
//cout << "Updated C_mat[" << r << "][" << r << "] = " << C_mat[r][r] << endl;

    // Edge Confining and Ionic Potential Terms:-->
    if(Parameters_BHZ_.canonical_==false){
        int rx, ry;
        ry=(r-spin*half_size_-orb)%(2*ly_);
        rx=(r-spin*half_size_-orb)/(2*ly_);

        if(Parameters_BHZ_.numberofSoftedges_==1){
            if (rx > lx_ - wx_){
                C_mat[r][r] += (Vo_ * 1.0 * (rx + wx_ - lx_) / (1.0 * (wx_ - 1))) * One_Complex;
                
            }
            if(rx == lx_- wx_){
                C_mat[r][r] += Parameters_BHZ_.helical_edge_pot * One_Complex;
                if(spin==0){
                    C_mat[r][r] += Parameters_BHZ_.Hfield * One_Complex;
                }
                else{
                    C_mat[r][r] += -1.0 * Parameters_BHZ_.Hfield * One_Complex;
                }
                    
            }
        }
        else{
            if (rx < wx_ - 1){
                C_mat[r][r] += (Vo_ * 1.0 * (wx_ - 1 - rx) / (1.0 * (wx_ - 1))) * One_Complex;
            }
            if (rx > lx_ - wx_){
                C_mat[r][r] += (Vo_ * 1.0 * (rx + wx_ - lx_) / (1.0 * (wx_ - 1))) * One_Complex;
            }
            if(rx == lx_- wx_ || rx == wx_-1){
                C_mat[r][r] += Parameters_BHZ_.helical_edge_pot * One_Complex;
            }
        }
    }

/*    int rx=(r-spin*half_size_-orb)/(2*ly_);
    cout<<"r_x="<<rx<<"     "<<C_mat[r][r]<<endl;
*/

    // Adding Onsite Mass Term:-->
	C_mat[r][r] += ( 1.0 * (C_ - 4.0 * D_) + 1.0 * (pow(-1.0, 1.0 * orb)) * (M_ - 4.0 * B_) ) * One_Complex;
}

void Hamiltonian_BHZ::addHoppingXTerm(int r, int rpx, int orb){
    // Adding NN Hopping connections:
    C_mat[r][rpx] = ( 1.0 * D_ + 1.0 * (pow(-1.0, 1.0 * orb)) * B_ )* One_Complex;
    C_mat[rpx][r] = conj(C_mat[r][rpx]);
}

void Hamiltonian_BHZ::addHoppingYTerm(int r, int rpy, int orb){
    // Adding NN Hopping connections:
    C_mat[r][rpy] = ( 1.0 * D_ + 1.0 * (pow(-1.0, 1.0 * orb)) * B_ ) * One_Complex;
    C_mat[rpy][r] = conj(C_mat[r][rpy]);
}

void Hamiltonian_BHZ::addTBC_HoppingYTerm(int r, int rpy, int orb, int spin){
    //Flux-> spin dependent: 
    //from (rx,Ly-1) -> (rx,0) multiplies e^{i*spin*Phi_y}
    C_mat[r][rpy] = ( 1.0 * (pow(-1.0, 1.0 * orb)) * B_ * One_Complex ) * exp( (pow(-1.0, 1.0*(1-spin))) * Iota_Complex * Parameters_BHZ_.Phi_y );
    C_mat[rpy][r] = conj(C_mat[r][rpy]);

    //Flux-> charge independent:
    //from (rx,Ly-1) -> (rx,0) multiplies e^{i*Phi_y}
//    C_mat[r][rpy] = ( 1.0 * (pow(-1.0, 1.0 * orb)) * B_ * One_Complex ) * exp( (-1.0)* Iota_Complex * Parameters_BHZ_.Phi_y );
//    C_mat[rpy][r] = conj(C_mat[r][rpy]);
}

void Hamiltonian_BHZ::addOrbMixingXTerm(int r, int rpx, int orb, int spin){
    // Adding NN Orbital-Mixing connections:
    if(spin == 0){
        C_mat[rpx][r] = -(1.0 * A_ / 2.0) * Iota_Complex;
        C_mat[r][rpx] = conj(C_mat[rpx][r]);
    }
    if(spin == 1){
        C_mat[rpx][r] = (1.0 * A_ / 2.0) * Iota_Complex;
        C_mat[r][rpx] = conj(C_mat[rpx][r]);
    }
}

void Hamiltonian_BHZ::addOrbMixingYTerm(int r, int rpy, int orb, int spin){
    // Adding NN Orbital-Mixing connections:
    if(orb == 0){
        C_mat[rpy][r] = (1.0 * A_ / 2.0) * One_Complex;
        C_mat[r][rpy] = conj(C_mat[rpy][r]);
    }
    if(orb == 1){
        C_mat[rpy][r] = -(1.0 * A_ / 2.0) * One_Complex;
        C_mat[r][rpy] = conj(C_mat[rpy][r]);        
    }
}

void Hamiltonian_BHZ::addTBC_OrbMixingYTerm(int r, int rpy, int orb, int spin){
    //Flux-> spin dependent:
    //from (rx,Ly-1) -> (rx,0) multiplies e^{i*spin*Phi_y}
    if(orb==0){
        C_mat[r][rpy] = ( (1.0 * A_ / 2.0) * One_Complex ) * exp( (pow(-1.0, 1.0*(1-spin))) * Iota_Complex * Parameters_BHZ_.Phi_y);
        C_mat[rpy][r] = conj(C_mat[r][rpy]);
    }
    if(orb==1){
        C_mat[r][rpy] = ( -(1.0 * A_ / 2.0) * One_Complex ) * exp( (pow(-1.0, 1.0*(1-spin))) * Iota_Complex * Parameters_BHZ_.Phi_y);
        C_mat[rpy][r] = conj(C_mat[r][rpy]);
    }

    //Flux-> charge dependent:
    //from (rx,Ly-1) -> (rx,0) multiplies e^{i*Phi_y}
/*    if(orb==0){
        C_mat[r][rpy] = ( (1.0 * A_ / 2.0) * One_Complex ) * exp( (-1.0) * Iota_Complex * Parameters_BHZ_.Phi_y);
        C_mat[rpy][r] = conj(C_mat[r][rpy]);
    }
    if(orb==1){
        C_mat[r][rpy] = ( -(1.0 * A_ / 2.0) * One_Complex ) * exp( (-1.0) * Iota_Complex * Parameters_BHZ_.Phi_y);
        C_mat[rpy][r] = conj(C_mat[r][rpy]);
    }
*/
}

void Hamiltonian_BHZ::addPinningTerm(int r, int orb, int spin){
    if(abs(Parameters_BHZ_.Pinfield) > 1e-3){
        C_mat[r][r] += pow(-1.0,1.0*spin)*Parameters_BHZ_.Pinfield;
    }
}

void Hamiltonian_BHZ::connectionMatrix(){

    for(int i=0;i<H_size;i++){
        for(int j=0;j<H_size;j++){
            C_mat[i][j]=Zero_Complex;
        }
    }

    int r, rpx, rpy, bar_rpx, bar_rpy;
    for (int spin = 0; spin < N_spin_; spin++){
        for(int orb = 0; orb <N_orbs_;orb++){
            for(int rx=0;rx<lx_;rx++){
                for(int ry=0;ry<ly_;ry++){

                    r   = makeIndex(rx, ry, orb, spin);
                    addOnsiteTerms(r, orb, spin);

                //    if(rx==lx_-1 && ry==ly_-1){
                //        addPinningTerm(r, orb, spin);
                //    }

                    if(rx != lx_-1){
                        rpx = makeIndex(rx+1, ry, orb, spin);
                        addHoppingXTerm(r, rpx, orb);

                        bar_rpx = makeIndex(rx+1, ry, bar_function(orb), spin);
                        addOrbMixingXTerm(r, bar_rpx, orb, spin);
                    }
                    if(rx == lx_-1 && Parameters_BHZ_.PBC_X==true){
                        assert(lx_ > 2 && "lx_ must be greater than 2 for PBC along-x.");
                        rpx = makeIndex(0, ry, orb, spin);
                        addHoppingXTerm(r, rpx, orb);

                        bar_rpx = makeIndex(0, ry, bar_function(orb), spin);
                        addOrbMixingXTerm(r, bar_rpx, orb, spin);
                    }
                    if(ry != ly_-1){
                        rpy = makeIndex(rx, ry+1, orb, spin);
                        addHoppingYTerm(r, rpy, orb);

                        bar_rpy = makeIndex(rx, ry+1, bar_function(orb), spin);
                        addOrbMixingYTerm(r, bar_rpy, orb, spin);
                    }
                    if(ry == ly_-1 && Parameters_BHZ_.PBC_Y==true){
                        assert(ly_ > 2 && "ly_ must be greater than 2 for PBC along-y.");

                        if(Parameters_BHZ_.TBC_Y==true){
                            rpy = makeIndex(rx, 0, orb, spin);
                            addTBC_HoppingYTerm(r, rpy, orb, spin);

                            bar_rpy = makeIndex(rx, 0, bar_function(orb), spin);
                            addTBC_OrbMixingYTerm(r, bar_rpy, orb, spin);
                        }
                        else{
                            rpy = makeIndex(rx, 0, orb, spin);
                            addHoppingYTerm(r, rpy, orb);

                            bar_rpy = makeIndex(rx, 0, bar_function(orb), spin);
                            addOrbMixingYTerm(r, bar_rpy, orb, spin);
                        }
                        
                    }
                }
            }
        }
    }

    cout<<"Hamiltonian matrix constructed succesfully"<<endl;

/*
    for(int i=0;i<H_size;i++){
        for(int j=0;j<H_size;j++){
            cout<<C_mat[i][j]<<" ";
        }
        cout<<endl;
    }
*/
}

void Hamiltonian_BHZ::Diagonalizer(){

    cout<<"Starting the diagonalizer"<<endl;
    std::vector<std::complex<double>> Ham_(H_size*H_size);
    //std::vector<lapack_complex_double> Ham_(H_size*H_size);

    
    //#pragma omp parallel for default(shared) 
    for (int i = 0; i < H_size; ++i) {
        for (int j = 0; j < H_size; ++j) {
            //Ham_[i*H_size + j] = lapack_make_complex_double(C_mat[i][j].real(), C_mat[i][j].imag());
            Ham_[j*H_size + i] = C_mat[i][j];
        }
    }

    cout<<"C_mat copied in Ham"<<endl;

    // LAPACK routine variables
    char jobz = 'V'; // Computing both eigenvalues and eigenvectors
    char uplo = 'L'; // Using the lower triangular part of the matrix
    int n = H_size;
    int lda = H_size;
    int info;

    std::vector<double> eigs_(H_size);

    std::vector<std::complex<double>> work(1);
    //std::vector<lapack_complex_double> work(1);
    std::vector<double> rwork(3 * n - 2);
    int lwork = -1;

    //querying with lwork=-1
    zheev_(&jobz, &uplo, &n, Ham_.data(), &lda, eigs_.data(), work.data(), &lwork, rwork.data(), &info);

    cout<<"Diagonalization begins"<<endl;
    // Set optimal workspace size
    lwork = static_cast<int>(work[0].real());
    //lwork = static_cast<int>(lapack_complex_double_real(work[0]));
    work.resize(lwork);

    // Perform the eigenvalue decomposition
    zheev_(&jobz, &uplo, &n, Ham_.data(), &lda, eigs_.data(), work.data(), &lwork, rwork.data(), &info);

    // Check for successful execution
    if (info != 0) {
        std::cerr << "LAPACK zheev_ failed with info=" << info << std::endl;
        return;
    }
    else{
        cout<<"Ham diagonalized succesfully"<<endl;
    }
    cout<<"Diagonalization finished"<<endl;
    for(int i=0;i<H_size;i++){
        evals_[i] = eigs_[i];
    }
    eigs_.clear();

    for(int i=0;i<H_size;i++){
        for(int j=0;j<H_size;j++){
            evecs_[i][j] = Ham_[j*H_size+i];
        }
    }
    Ham_.clear();

/*
    for(int i=0;i<H_size;i++){
        for(int j=0;j<H_size;j++){
            evecs_[i * H_size + j] = Ham_[i*H_size+j];
        }
    }
    Ham_.clear();
*/

    string Evals_out="Eigenvalues.txt";
    ofstream Evals_file_out(Evals_out.c_str());

    for(int n=0;n<evals_.size();n++){
        Evals_file_out<<n<<"    "<<evals_[n]<<endl;
    }

    if(Parameters_BHZ_.save_evecs==true){
        string file_evecs="Eigenvectors.txt";
        ofstream file_evecs_out(file_evecs.c_str());

        for(int i=0;i<H_size;i++){
            for(int j=0;j<H_size;j++){
                file_evecs_out<<evecs_[i][j]<<" ";
            }
            file_evecs_out<<endl;
        }
    }
}

void Hamiltonian_BHZ::updownOccupiedSubspace(){

    int N_occupied;
    if(Parameters_BHZ_.canonical_){
        N_occupied = Parameters_BHZ_.total_particles_;
    }
    else{
        cerr<<"N_occupied not yet defined for GCE"<<endl;
    }

    Mat_1_int up_occ_indices,dn_occ_indices;
    up_occ_indices.clear();     dn_occ_indices.clear();

    for(int n=0;n<N_occupied;n++){
        double proj_up = 0.0;
        double proj_dn = 0.0;

        for(int r=0;r<half_size_;r++){//spin-up Connection matrix block
            proj_up += norm(evecs_[r][n]);
        }
        for(int r=half_size_;r<H_size;r++){//spin-dn Connection matrix block
            proj_dn += norm(evecs_[r][n]);
        }

        if(proj_up >= proj_dn){
            up_occ_indices.push_back(n);
        }
        else{
            dn_occ_indices.push_back(n);
        }
    }

    //Constructing matrices for occupied subspace of spin-up and spin-down eigenvectors (row-major):
    Mat_2_Complex_doub spin_up_HS,spin_dn_HS;
    spin_up_HS.resize(up_occ_indices.size());
    spin_dn_HS.resize(dn_occ_indices.size());

    for(int i=0;i<up_occ_indices.size();i++){
        int band_index = up_occ_indices[i];
        spin_up_HS[i].resize(H_size);
        for(int r=0;r<H_size;r++){
            spin_up_HS[i][r] = evecs_[r][band_index];
        }
    }

    for(int i=0;i<dn_occ_indices.size();i++){
        int band_index = dn_occ_indices[i];
        spin_dn_HS[i].resize(H_size);
        for(int r=0;r<H_size;r++){
            spin_dn_HS[i][r] = evecs_[r][band_index];
        }
    }

    string file_evecs_up="Occupied_evecs_spin_up.txt";
    ofstream file_evecs_up_out(file_evecs_up.c_str());

    string file_evecs_dn="Occupied_evecs_spin_dn.txt";
    ofstream file_evecs_dn_out(file_evecs_dn.c_str());

    for(int i=0;i<up_occ_indices.size();i++){
        for(int r=0;r<H_size;r++){
            file_evecs_up_out<<spin_up_HS[i][r]<<" ";
        }
        file_evecs_up_out<<endl;
    }

    for(int i=0;i<dn_occ_indices.size();i++){
        for(int r=0;r<H_size;r++){
            file_evecs_dn_out<<spin_dn_HS[i][r]<<" ";
        }
        file_evecs_dn_out<<endl;
    }

    cout<<"total # of occupied spin-up evecs ="<<up_occ_indices.size()<<endl;
    cout<<"total # of occupied spin-dn evecs ="<<dn_occ_indices.size()<<endl;
}

#endif
