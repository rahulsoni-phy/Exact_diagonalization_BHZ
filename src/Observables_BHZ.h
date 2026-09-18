#include "tensors.h"
#include "Parameters_BHZ.h"
#include "Hamiltonian_BHZ.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef Observables_BHZ_class
#define Observables_BHZ_class
#define PI acos(-1.0)

class Observables_BHZ{

    public:
    Observables_BHZ(Parameters_BHZ &Parameters_BHZ__, Hamiltonian_BHZ &Hamiltonian_BHZ__): 
        Parameters_BHZ_(Parameters_BHZ__),Hamiltonian_BHZ_(Hamiltonian_BHZ__){
            Initialize();
    }

    Parameters_BHZ &Parameters_BHZ_;
    Hamiltonian_BHZ &Hamiltonian_BHZ_;

    void Initialize();
    double fermifunction(double en_, double mu_);
    double chemicalpotential(int particles_);
    void calculateDOS();
    complex<double> calculateAmplitude(int state1, int pos1, int state2, int pos2);
    void calculateLDOE();
    void calculateSpinCurrents();
//    void calculateBmat();
    void calculateAkxw();
    void calculateAkyw();
    void calculateLDOS();
//    void calculateYMomentumResolvedAkxw();
    void calculateWaveFunctions();
    void calculateAkywNew();
    void calculateEdgeConductance();

    int lx_, ly_, H_size, MHS;
    int N_cells_,N_orbs_, N_spin_;
    double mu;
    double one_by_PI_=1.0/(1.0*PI);
    double w_min,w_max,dw,eta;
    int w_size;

};

void Observables_BHZ::Initialize(){

    lx_ = Parameters_BHZ_.Lx_;
    ly_ = Parameters_BHZ_.Ly_;
    
    N_cells_ = Parameters_BHZ_.Total_Cells;
    N_orbs_ = Parameters_BHZ_.N_Orbs;
    N_spin_ = 2;
    H_size = Parameters_BHZ_.Ham_Size;
    MHS = Hamiltonian_BHZ_.half_size_;

    w_min = Parameters_BHZ_.w_min;
    w_max = Parameters_BHZ_.w_max;

    dw = Parameters_BHZ_.dw_;
    eta = Parameters_BHZ_.eta_;

    if(Parameters_BHZ_.canonical_){
        mu=chemicalpotential(Parameters_BHZ_.total_particles_);
        cout<<"Canonically calculated mu = "<<mu<<endl;
    }
    else{
        mu=Parameters_BHZ_.mu_fixed;
        cout<<"Grand-canonical fixed mu = "<<mu<<endl;
    }

    w_size = (int) ( (w_max - w_min)/dw );
}

double Observables_BHZ::chemicalpotential(int particles_){
    double mu_temp, eps_, dmu_by_dN, N_temp, dmu_by_dN_min, Ne_;
    eps_=1e-4;
    mu_temp=Hamiltonian_BHZ_.evals_[0];
    N_temp=100000;
    Ne_=1.0*particles_;

    int iters=0;

    dmu_by_dN = 0.01*( Hamiltonian_BHZ_.evals_[H_size-1] - Hamiltonian_BHZ_.evals_[0] )*( 1.0/(1.0*H_size) );
    dmu_by_dN_min = 0.0001*( Hamiltonian_BHZ_.evals_[H_size-1] - Hamiltonian_BHZ_.evals_[0] )*( 1.0/(1.0*H_size) );

    while( abs(N_temp - Ne_) > eps_){
        N_temp=0;
        for(int n=0;n<H_size;n++){
            N_temp += fermifunction(Hamiltonian_BHZ_.evals_[n],mu_temp);
        }

        mu_temp = mu_temp + dmu_by_dN*(Ne_-N_temp);
        iters++;

        if(iters%1000==0){
            dmu_by_dN = (1000.0/(10.0*iters))*dmu_by_dN;
            dmu_by_dN = max(dmu_by_dN, dmu_by_dN_min);
        }
    }

    cout<<"Calculated # of particles = "<<N_temp<<endl;

    return mu_temp;
}

double Observables_BHZ::fermifunction(double en_, double mu_){
	double ffn, temp;
	temp = Parameters_BHZ_.temp_;
	ffn = ((1.0) / (1.0 + exp((en_ - mu_) / temp)));
	return ffn;
}

void Observables_BHZ::calculateDOS(){
    double value, w;

    string DOS_out="Density_of_states.txt";
    ofstream DOS_file_out(DOS_out.c_str());

    w=w_min;
    while(w<=w_max){
        value=0.0;
        DOS_file_out<< w<<"     ";

        for(int n=0;n<H_size;n++){
            value=value+(1.0/(4.0*lx_*ly_))*(one_by_PI_)*((eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)));
        }
        DOS_file_out<<value<<endl;
        w=w+dw;
    }
}

complex<double> Observables_BHZ::calculateAmplitude(int state1, int pos1, int state2, int pos2){
    complex<double> value;
    value = (conj(Hamiltonian_BHZ_.evecs_[pos1][state1]))*(Hamiltonian_BHZ_.evecs_[pos2][state2]);
    return value;
}

void Observables_BHZ::calculateLDOE(){

    Mat_4_Complex_doub ChargeDensity;
    ChargeDensity.resize(lx_);
    for(int rx=0;rx<lx_;rx++){
        ChargeDensity[rx].resize(ly_);
        for(int ry=0;ry<ly_;ry++){
            ChargeDensity[rx][ry].resize(N_orbs_);
            for(int orb=0;orb<N_orbs_;orb++){
                ChargeDensity[rx][ry][orb].resize(N_spin_);
                for(int spin=0;spin<N_spin_;spin++){
                    ChargeDensity[rx][ry][orb][spin]=Zero_Complex;
                }
            }
        }
    }

    complex<double> Tot_num_ples, local_charge_density;
    Tot_num_ples=Zero_Complex;

    int r;
    for(int rx=0;rx<lx_;rx++){
        for(int ry=0;ry<ly_;ry++){
            for(int orb=0;orb<N_orbs_;orb++){
                for(int spin=0;spin<N_spin_;spin++){
                    r = Hamiltonian_BHZ_.makeIndex(rx, ry, orb, spin);

                    local_charge_density = Zero_Complex;
                    for(int n=0;n<H_size;n++){
                        local_charge_density += calculateAmplitude(n, r, n, r)*fermifunction(Hamiltonian_BHZ_.evals_[n],mu);
                    }
                    ChargeDensity[rx][ry][orb][spin] = local_charge_density;
                    Tot_num_ples += ChargeDensity[rx][ry][orb][spin];
                }
            }
        }
    }
    cout<<"Total # of particles = "<<Tot_num_ples.real()<<endl;

    double inv_ly=(1.0/(1.0*ly_));

    string file_Avg_LDOE="LDOE_plus_sz_along_x.txt";
    ofstream file_Avg_CD_out(file_Avg_LDOE.c_str());
    if (!file_Avg_CD_out.is_open()){
        cerr << "Error: Could not open file " << file_Avg_LDOE << endl;
        return;
    }

    file_Avg_CD_out <<"#rx       LDOE_Orb-0    LDOE_Orb-1    Sz_Orb-0    Sz_Orb-1"<< endl;

    complex<double> CD_s_up,CD_s_dn,CD_p_up,CD_p_dn;
    for(int rx=0;rx<lx_;rx++){
        CD_s_up = Zero_Complex;
        CD_s_dn = Zero_Complex;
        CD_p_up = Zero_Complex;
        CD_p_dn = Zero_Complex;

        for(int ry=0;ry<ly_;ry++){
            CD_s_up += inv_ly*ChargeDensity[rx][ry][0][0];
            CD_s_dn += inv_ly*ChargeDensity[rx][ry][0][1];
            CD_p_up += inv_ly*ChargeDensity[rx][ry][1][0];
            CD_p_dn += inv_ly*ChargeDensity[rx][ry][1][1];
        }
        file_Avg_CD_out<<rx<<"  "<<CD_s_up.real()+CD_s_dn.real()<<"    "<<CD_p_up.real()+CD_p_dn.real()<<"      "<<
                                    (CD_s_up.real()-CD_s_dn.real())*0.5<<"       "<<0.5*(CD_p_up.real()-CD_p_dn.real())<<endl;
    }

    if(!Parameters_BHZ_.canonical_){
    string file_Bulk_Edge_LDOE="Bulk_edge_particles.txt";
    ofstream file_BE_out(file_Bulk_Edge_LDOE.c_str());
    if (!file_BE_out.is_open()) {
        cerr << "Error: Could not open file " << file_Bulk_Edge_LDOE << endl;
        return;
    }

    complex<double> Tot_edge_ples_orb_s,Tot_edge_ples_orb_p;
    Tot_edge_ples_orb_s=Zero_Complex;
    Tot_edge_ples_orb_p=Zero_Complex;

    complex<double> Tot_edge_ples, Tot_bulk_ples;
    Tot_edge_ples=Zero_Complex;
    Tot_bulk_ples=Zero_Complex;

    for(int rx=0;rx<lx_;rx++){
        for(int ry=0;ry<ly_;ry++){
            if(Parameters_BHZ_.numberofSoftedges_==1){
                if(rx>lx_-Hamiltonian_BHZ_.wx_){
                    Tot_edge_ples_orb_s += ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1];
                    Tot_edge_ples_orb_p += ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                    Tot_edge_ples=Tot_edge_ples_orb_s + Tot_edge_ples_orb_p;
                }
                else{
                    Tot_bulk_ples+=ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1] + ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                }
            }
            else{
                if(rx<Hamiltonian_BHZ_.wx_-1 || rx>lx_-Hamiltonian_BHZ_.wx_){
                    Tot_edge_ples_orb_s += ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1];
                    Tot_edge_ples_orb_p += ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                    Tot_edge_ples=Tot_edge_ples_orb_s + Tot_edge_ples_orb_p;
                }
                else{
                    Tot_bulk_ples+=ChargeDensity[rx][ry][0][0] + ChargeDensity[rx][ry][0][1] + ChargeDensity[rx][ry][1][0] + ChargeDensity[rx][ry][1][1];
                }
            }
        }
    }

    file_BE_out<<"Edge_ples_orb_s= "<<Tot_edge_ples_orb_s.real()<<endl;
    file_BE_out<<"Edge_ples_orb_p= "<<Tot_edge_ples_orb_p.real()<<endl;
    file_BE_out<<"Total_edge_ples= "<<Tot_edge_ples.real()<<endl;
    file_BE_out<<"Total_bulk_ples= "<<Tot_bulk_ples.real()<<endl;
    file_BE_out<<"Total_num_ples= "<<Tot_num_ples.real()<<endl;
    }

    if(Parameters_BHZ_.TBC_Y){
        Hamiltonian_BHZ_.updownOccupiedSubspace();
        
        string file_flux="Flux_avgd_observes.txt";
        ofstream file_flux_out(file_flux.c_str());

        complex<double> halfocc_s_up, halfocc_s_dn, halfocc_p_up, halfocc_p_dn;
        complex<double> halfSz_s, halfSz_p;

        int Lxby2;
        if(lx_%2==0){Lxby2 = lx_/2;}
        else{Lxby2 = (lx_-1)/2;}

        halfocc_s_up = Zero_Complex;    halfocc_s_dn = Zero_Complex;
        halfocc_p_up = Zero_Complex;    halfocc_p_dn = Zero_Complex;
        halfSz_s = Zero_Complex;         halfSz_p = Zero_Complex;
        for(int rx=0; rx<Lxby2; rx++){            
            for(int ry=0; ry< ly_; ry++){
                halfocc_s_up += ChargeDensity[rx][ry][0][0];
                halfocc_s_dn += ChargeDensity[rx][ry][0][1];
                halfocc_p_up += ChargeDensity[rx][ry][1][0];
                halfocc_p_dn += ChargeDensity[rx][ry][1][1];
            }
        }
        halfSz_s = 0.5*(halfocc_s_up - halfocc_s_dn);
        halfSz_p = 0.5*(halfocc_p_up - halfocc_p_dn);

        file_flux_out<<"Total_ns_up_half_cylinder= "<<halfocc_s_up.real()<<endl;
        file_flux_out<<"Total_ns_dn_half_cylinder= "<<halfocc_s_dn.real()<<endl;
        file_flux_out<<"Total_np_up_half_cylinder= "<<halfocc_p_up.real()<<endl;
        file_flux_out<<"Total_np_dn_half_cylinder= "<<halfocc_p_dn.real()<<endl;

        file_flux_out<<"Total_Sz_s_half_cylinder= "<<halfSz_s.real()<<endl;
        file_flux_out<<"Total_Sz_p_half_cylinder= "<<halfSz_p.real()<<endl;

        file_flux_out<<"Total_Sz_half_cylinder= "<<(halfSz_s.real() + halfSz_p.real())<<endl;

        int i0=H_size/2-1;   
        double E_single_gap = Hamiltonian_BHZ_.evals_[i0+1] -Hamiltonian_BHZ_.evals_[i0];
        file_flux_out<<"Single_particle_Gap= "<<E_single_gap<<endl;

        int qmax = 3;
        int numoftargetstate=2*qmax + 1;

        Mat_1_int band_index;
        band_index.reserve(numoftargetstate);
        for (int q=-qmax;q<=qmax;q++){      band_index.push_back(i0 + q);   }

        Mat_2_doub amplitude_vec;
        amplitude_vec.resize(lx_);
        for(int rx=0;rx<lx_;rx++){
            amplitude_vec[rx].resize(numoftargetstate);
        }
        for(int k=0;k<numoftargetstate;k++){
            for(int rx=0;rx<lx_;rx++){
                amplitude_vec[rx][k] = 0.0;
                for(int ry=0;ry<ly_;ry++){
                    for(int orb=0;orb<N_orbs_;orb++){
                        for(int spin=0;spin<N_spin_;spin++){
                            r = Hamiltonian_BHZ_.makeIndex(rx, ry, orb, spin);

                            amplitude_vec[rx][k] += inv_ly*norm(Hamiltonian_BHZ_.evecs_[r][band_index[k]]);
                        }
                    }
                }
            }
        }

        string file_wavefn="Wavefunction_amplitudes_for_q_m2_to_p2.txt";
        ofstream file_out_wavefn(file_wavefn.c_str());
        file_out_wavefn<<"#rx   A(i0-3)     A(i0-2)     A(i0-1)     A(i0)       A(i0+1)     A(i0+2)     A(i0+3)"<<endl;
        for(int rx=0;rx<lx_;rx++){
            file_out_wavefn<<rx<<"  ";
            for(int k=0;k<numoftargetstate;k++){
                file_out_wavefn<<amplitude_vec[rx][k]<<"        ";
            }
            file_out_wavefn<<endl;
        }

    }

}


void Observables_BHZ::calculateSpinCurrents(){
    Mat_4_Complex_doub SpinCurrent_UP, SpinCurrent_DN;

    SpinCurrent_UP.resize(N_cells_);
    SpinCurrent_DN.resize(N_cells_);
    for(int cell1=0;cell1<N_cells_;cell1++){
        SpinCurrent_UP[cell1].resize(N_cells_);
        SpinCurrent_DN[cell1].resize(N_cells_);

        for(int cell2=0;cell2<N_cells_;cell2++){
            SpinCurrent_UP[cell1][cell2].resize(N_orbs_);
            SpinCurrent_DN[cell1][cell2].resize(N_orbs_);

            for(int orb1=0;orb1<N_orbs_;orb1++){
                SpinCurrent_UP[cell1][cell2][orb1].resize(N_orbs_);
                SpinCurrent_DN[cell1][cell2][orb1].resize(N_orbs_);

                for(int orb2=0;orb2<N_orbs_;orb2++){
                    SpinCurrent_UP[cell1][cell2][orb1][orb2]=Zero_Complex;
                    SpinCurrent_DN[cell1][cell2][orb1][orb2]=Zero_Complex;
                }
            }
        }
    }

    int r1, r2, cell1, cell2;
//    complex<double> local_spinup_curr, local_spindn_curr;
//    #pragma omp parallel for collapse(6) private(r1, r2, cell1, cell2)
    for(int r1x=0;r1x<lx_;r1x++){
        for(int r1y=0;r1y<ly_;r1y++){
            cell1 = r1y + ly_*r1x;

            for(int orb1=0;orb1<N_orbs_;orb1++){
                r1 = Hamiltonian_BHZ_.makeIndex(r1x, r1y, orb1, 0);

                for(int r2x=0;r2x<lx_;r2x++){
                    for(int r2y=0;r2y<ly_;r2y++){
                        cell2 = r2y + ly_*r2x;

                        for(int orb2=0;orb2<N_orbs_;orb2++){
                            r2 = Hamiltonian_BHZ_.makeIndex(r2x, r2y, orb2, 0);

                            if(cell1 < cell2){
                                if(Hamiltonian_BHZ_.C_mat[r1][r2].real()!=0 || Hamiltonian_BHZ_.C_mat[r1][r2].imag()!=0){
                                    for(int n=0;n<H_size;n++){
                                        SpinCurrent_UP[cell1][cell2][orb1][orb2] += 
                                        0.5*Iota_Complex*( ( Hamiltonian_BHZ_.C_mat[r2][r1]*calculateAmplitude(n, r2, n, r1) ) - 
                                        ( Hamiltonian_BHZ_.C_mat[r1][r2]*calculateAmplitude(n, r1, n, r2) ) )*(fermifunction(Hamiltonian_BHZ_.evals_[n],mu) );
                                    }
                                }
                                if(Hamiltonian_BHZ_.C_mat[r1+MHS][r2+MHS].real()!=0 || Hamiltonian_BHZ_.C_mat[r1+MHS][r2+MHS].imag()!=0){
                                    for(int n=0;n<H_size;n++){
                                        SpinCurrent_DN[cell1][cell2][orb1][orb2] 
                                        += 0.5*Iota_Complex*( ( Hamiltonian_BHZ_.C_mat[r2+MHS][r1+MHS]*calculateAmplitude(n, r2+MHS, n, r1+MHS) )- 
                                        ( Hamiltonian_BHZ_.C_mat[r1+MHS][r2+MHS]*calculateAmplitude(n, r1+MHS, n, r2+MHS) ) )*(fermifunction(Hamiltonian_BHZ_.evals_[n],mu) );
                                    }
                                }
                            //    cout<<cell1<<"  "<<cell2<<" "<<orb1<<"  "<<orb2<<"  "<<SpinCurrent_UP[cell1][cell2][orb1][orb2].real()<<"   "
                            //    <<SpinCurrent_DN[cell1][cell2][orb1][orb2].real()<<endl;
                            }
                        }
                    }
                }
            }
        }
    }

    string file_spin_current="Total_spin_current_at_each_cell_link.txt";
    ofstream file_spin_current_out(file_spin_current.c_str());
    if (!file_spin_current_out.is_open()) {
        cerr << "Error: Could not open file " << file_spin_current << endl;
        return;
    }

    for(int r1x=0;r1x<lx_;r1x++){
        for(int r1y=0;r1y<ly_;r1y++){
            cell1 = r1y + ly_*r1x;

            for(int r2x=0;r2x<lx_;r2x++){
                for(int r2y=0;r2y<ly_;r2y++){
                    cell2 = r2y + ly_*r2x;

                    if(cell1<cell2){
                        if(SpinCurrent_UP[cell1][cell2][0][0].real()!=0 || SpinCurrent_UP[cell1][cell2][1][1].real()!=0){
                            file_spin_current_out<<cell1<<" "<<cell2<<" "
                            <<SpinCurrent_UP[cell1][cell2][0][0].real()-SpinCurrent_DN[cell1][cell2][0][0].real()<<"  "
                            <<SpinCurrent_UP[cell1][cell2][1][1].real()-SpinCurrent_DN[cell1][cell2][1][1].real()<<endl;
                        }
                    }
                }
            }
        }
    }

    string file_avg_spin_current="Avg_spin_current_along_ry_vs_rx.txt";
    ofstream file_avg_spin_current_out(file_avg_spin_current.c_str());
    if (!file_avg_spin_current_out.is_open()) {
        cerr << "Error: Could not open file " << file_avg_spin_current << endl;
        return;
    }

    file_avg_spin_current_out<<"#rx     #(s->s),up      #(p->p),up      #(s->s),dn      #(p->p),dn"<<endl;
    
    complex<double> SC_ss_up,SC_pp_up,SC_ss_dn,SC_pp_dn;
    for(int r1x=0;r1x<lx_;r1x++){
        SC_ss_up=Zero_Complex;  SC_pp_up=Zero_Complex;
        SC_ss_dn=Zero_Complex;  SC_pp_dn=Zero_Complex;

        for(int r1y=0;r1y<ly_;r1y++){
            cell1 = r1y + ly_*r1x;

            for(int r2y=0;r2y<ly_;r2y++){
                cell2 = r2y + ly_*r1x;

                if(cell1<cell2){
                    SC_ss_up += SpinCurrent_UP[cell1][cell2][0][0];
                    SC_pp_up += SpinCurrent_UP[cell1][cell2][1][1];

                    SC_ss_dn += SpinCurrent_DN[cell1][cell2][0][0];
                    SC_pp_dn += SpinCurrent_DN[cell1][cell2][1][1];
                }
            }
        }
        file_avg_spin_current_out<<r1x<<"   "<<SC_ss_up.real()<<"     "<<SC_pp_up.real()<<"     "<<SC_ss_dn.real()<<"   "<<SC_pp_dn.real()<<endl;
    }
    
}

void Observables_BHZ::calculateLDOS(){

    string file_ldos("LDOS_along_x.txt");
    ofstream file_ldos_out(file_ldos.c_str());

    double w=0.0;
    int r0_up, r0_dn, r1_up, r1_dn;
    complex<double> ldos_0_up, ldos_0_dn, ldos_1_up, ldos_1_dn;

    for(int rx=0;rx<lx_;rx++){
        for(int om=0;om<w_size;om++){
            w=w_min+om*dw;
            ldos_0_up=Zero_Complex;
            ldos_0_dn=Zero_Complex;
            ldos_1_up=Zero_Complex;
            ldos_1_dn=Zero_Complex;

            for(int ry=0;ry<ly_;ry++){
                r0_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 0);
                r0_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 1);
                r1_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 0);
                r1_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 1);

                for(int n=0;n<H_size;n++){
                    ldos_0_up += (one_by_PI_)*(calculateAmplitude(n, r0_up, n, r0_up))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                    ldos_0_dn += (one_by_PI_)*(calculateAmplitude(n, r0_dn, n, r0_dn))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                    ldos_1_up += (one_by_PI_)*(calculateAmplitude(n, r1_up, n, r1_up))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                    ldos_1_dn += (one_by_PI_)*(calculateAmplitude(n, r1_dn, n, r1_dn))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                }
                    
            }
            file_ldos_out<<rx<<"    "<<w<<"     "<<ldos_0_up.real()<<"   "<<ldos_0_dn.real()<<"     "<<ldos_1_up.real()<<"  "<<ldos_1_dn.real()<<endl;
        }
        file_ldos_out<<endl;
    }

}

void Observables_BHZ::calculateAkxw(){

    Mat_3_Complex_doub B_mat;
    B_mat.resize(lx_);
    for(int r1x=0;r1x<lx_;r1x++){
        B_mat[r1x].resize(lx_);
        for(int r2x=0;r2x<lx_;r2x++){
            B_mat[r1x][r2x].resize(w_size);
            for(int om=0;om<w_size;om++){
                B_mat[r1x][r2x][om]=Zero_Complex;
            }
        }
    }

    double w=0.0;
    int r1,r2;

    for(int r1x=0;r1x<lx_;r1x++){
        for(int r2x=0;r2x<lx_;r2x++){
            for(int om=0;om<w_size;om++){
                w=w_min+om*dw;

                for(int ry=0;ry<ly_;ry++){
                    for(int orb=0;orb<N_orbs_;orb++){
                        for(int spin=0;spin<N_spin_;spin++){
                            r1 = Hamiltonian_BHZ_.makeIndex(r1x, ry, orb, spin);
                            r2 = Hamiltonian_BHZ_.makeIndex(r2x, ry, orb, spin);

                            for(int n=0;n<H_size;n++){
                                    B_mat[r1x][r2x][om] += (one_by_PI_)*(calculateAmplitude(n, r1, n, r2))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                            }

                        }
                    }
                }
            }
        }
    }
    

    string file_Akxw("Akxw.txt");
    ofstream file_Akxw_out(file_Akxw.c_str());
    
    double kx,ky, momentum_step;
    complex<double> Akxw_;

    int kx_min,kx_max,kx_ind;
    if(Parameters_BHZ_.PBC_X){
        kx_min=-lx_/2;       kx_max=lx_/2;
        momentum_step = ((2.0*PI)/(1.0*lx_));
    }
    else{
        kx_min=1;       kx_max=lx_;
        momentum_step = ((1.0*PI)/(1.0*lx_+1.0));
    }

    for(kx_ind=kx_min;kx_ind<=kx_max;kx_ind++){
        kx=kx_ind*momentum_step;
        for(int om=0;om<w_size;om++){
            w=w_min+om*dw;
            Akxw_ = Zero_Complex;

            for(int r1x=0;r1x<lx_;r1x++){
                for(int r2x=0;r2x<lx_;r2x++){
                    for(int ry=0;ry<ly_;ry++){

                        if(Parameters_BHZ_.PBC_X){
                            Akxw_ += (1.0/(1.0*lx_))*( exp(Iota_Complex*kx*(1.0*(r1x-r2x)))*B_mat[r1x][r2x][om] );
                        }
                        else{
                            Akxw_ += (2.0/((lx_+1)*1.0))*( sin(kx*(1.0*r1x+1.0))*sin(kx*(1.0*r2x+1.0))*B_mat[r1x][r2x][om] );
                        }
                    }                    
                }
            }
            file_Akxw_out<<kx<<"    "<<w<<"         "<<Akxw_.real()<<endl;
        }
        file_Akxw_out<<endl;
    }

}

void Observables_BHZ::calculateAkyw(){
    Mat_4_Complex_doub B_s_mat,B_p_mat,B_mat;
    B_s_mat.resize(ly_);    B_p_mat.resize(ly_);
    B_mat.resize(ly_);

    for(int r1y=0;r1y<ly_;r1y++){
        B_s_mat[r1y].resize(ly_);    B_p_mat[r1y].resize(ly_);
        B_mat[r1y].resize(ly_);
        for(int r2y=0;r2y<ly_;r2y++){
            B_s_mat[r1y][r2y].resize(N_spin_);      B_p_mat[r1y][r2y].resize(N_spin_);
            B_mat[r1y][r2y].resize(N_spin_);
            for(int sp=0;sp<N_spin_;sp++){
                B_s_mat[r1y][r2y][sp].resize(w_size);       B_p_mat[r1y][r2y][sp].resize(w_size);
                B_mat[r1y][r2y][sp].resize(w_size);
                for(int om=0;om<w_size;om++){
                    B_s_mat[r1y][r2y][sp][om]=Zero_Complex;     B_p_mat[r1y][r2y][sp][om]=Zero_Complex;   
                    B_mat[r1y][r2y][sp][om]=Zero_Complex;
                }                
            }
        }
    }

    double w=0.0;
    int r1,r2;

    for(int r1y=0;r1y<ly_;r1y++){
        for(int r2y=0;r2y<ly_;r2y++){
            for(int spin=0;spin<N_spin_;spin++){
                for(int om=0;om<w_size;om++){
                    w=w_min+om*dw;

                    for(int orb=0;orb<N_orbs_;orb++){
                        for(int rx=0;rx<lx_;rx++){
                            
                            r1 = Hamiltonian_BHZ_.makeIndex(rx, r1y, orb, spin);
                            r2 = Hamiltonian_BHZ_.makeIndex(rx, r2y, orb, spin);

                            for(int n=0;n<H_size;n++){
                                    B_mat[r1y][r2y][spin][om] += (one_by_PI_)*(calculateAmplitude(n, r1, n, r2))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                            }

                            if(orb==0){
                                for(int n=0;n<H_size;n++){
                                    B_s_mat[r1y][r2y][spin][om] += (one_by_PI_)*(calculateAmplitude(n, r1, n, r2))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                                }
                            }

                            if(orb==1){
                                for(int n=0;n<H_size;n++){
                                    B_p_mat[r1y][r2y][spin][om] += (one_by_PI_)*(calculateAmplitude(n, r1, n, r2))*
                                        ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    string file_Akyw("Akyw.txt");
    ofstream file_Akyw_out(file_Akyw.c_str());
    file_Akyw_out<<"#ky"<<" "<<"omega"<<"   "<<"Akyw_up"<<"     "<<"Akyw_dn"<<"     "<<"Akyw_tot"<<endl;

    string file_Akyw_s_spin_resolved("Akyw_s_spin_resolved.txt");
    ofstream file_Akyw_s_sr_out(file_Akyw_s_spin_resolved.c_str());
    file_Akyw_s_sr_out<<"#ky"<<" "<<"omega"<<"   "<<"Akyw_s_up"<<"  "<<"Akyw_s_dn"<<"  "<<"Akyw_s_tot"<<endl;

    string file_Akyw_p_spin_resolved("Akyw_p_spin_resolved.txt");
    ofstream file_Akyw_p_sr_out(file_Akyw_p_spin_resolved.c_str());
    file_Akyw_p_sr_out<<"#ky"<<" "<<"omega"<<"   "<<"Akyw_p_up"<<"  "<<"Akyw_p_dn"<<"  "<<"Akyw_p_tot"<<endl;

    complex<double> Akyw_s,Akyw_s_up,Akyw_s_dn;
    complex<double> Akyw_p,Akyw_p_up,Akyw_p_dn;
    complex<double> Akyw_,Akyw_up,Akyw_dn;

    double momentum_step, ky;
    int ky_ind,ky_min,ky_max;   
    if(Parameters_BHZ_.PBC_Y){
        ky_min=0;       ky_max=ly_;
        momentum_step = ((2.0*PI)/(1.0*ly_));
    }
    else{
        ky_min=1;       ky_max=ly_;
        momentum_step = ((1.0*PI)/(1.0*ly_+1.0));
    }

    for(ky_ind=ky_min;ky_ind<=ky_max;ky_ind++){
        ky=ky_ind*momentum_step;

        for(int om=0;om<w_size;om++){
            w=w_min+om*dw;
            Akyw_s=Zero_Complex;Akyw_s_up=Zero_Complex;Akyw_s_dn=Zero_Complex;
            Akyw_p=Zero_Complex;Akyw_p_up=Zero_Complex;Akyw_p_dn=Zero_Complex;
            Akyw_=Zero_Complex;Akyw_up=Zero_Complex;Akyw_dn=Zero_Complex;

            for(int r1y=0;r1y<ly_;r1y++){
                for(int r2y=0;r2y<ly_;r2y++){

                        if(Parameters_BHZ_.PBC_Y){
                            Akyw_s_up += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_s_mat[r1y][r2y][0][om] );
                            Akyw_s_dn += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_s_mat[r1y][r2y][1][om] );
                            Akyw_s += (1.0/(2.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))* (B_s_mat[r1y][r2y][0][om] + B_s_mat[r1y][r2y][1][om]) );

                            Akyw_p_up += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_p_mat[r1y][r2y][0][om] );
                            Akyw_p_dn += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_p_mat[r1y][r2y][1][om] );
                            Akyw_p += (1.0/(2.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))* (B_p_mat[r1y][r2y][0][om] + B_p_mat[r1y][r2y][1][om]) );

                            Akyw_up += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_mat[r1y][r2y][0][om] );
                            Akyw_dn += (1.0/(1.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))*B_mat[r1y][r2y][1][om] );
                            Akyw_ += (1.0/(2.0*ly_*lx_))*( exp(Iota_Complex*ky*(1.0*(r1y-r2y)))* (B_mat[r1y][r2y][0][om] + B_mat[r1y][r2y][1][om]) );
                        }
                        else{
                            Akyw_ += (2.0/((ly_+1)*1.0))*( sin(ky*(1.0*r1y+1.0))*sin(ky*(1.0*r2y+1.0))* (B_mat[r1y][r2y][0][om] + B_mat[r1y][r2y][1][om]) );
                        }
                }
            }
            file_Akyw_out<<ky<<"    "<<w<<"         "<<Akyw_up.real()<<"    "<<Akyw_dn.real()<<"    "<<Akyw_.real()<<endl;
            file_Akyw_p_sr_out<<ky<<"    "<<w<<"         "<<Akyw_p_up.real()<<"     "<<Akyw_p_dn.real()<<"      "<<Akyw_p.real()<<endl;
            file_Akyw_s_sr_out<<ky<<"    "<<w<<"         "<<Akyw_s_up.real()<<"     "<<Akyw_s_dn.real()<<"      "<<Akyw_s.real()<<endl;
        }
        file_Akyw_s_sr_out<<endl;
        file_Akyw_p_sr_out<<endl;
        file_Akyw_out<<endl;
    }

}


void Observables_BHZ::calculateAkywNew(){
    //Evaluating the Lorentzian for each eigen-index and saving it in a 2D vector:
    Mat_2_doub lorentzian_;
    lorentzian_.resize(H_size);
    for(int n=0;n<H_size;n++){
        lorentzian_[n].resize(w_size);
        for(int om=0;om<w_size;om++){
            lorentzian_[n][om] = 0.0;
        }
    }

    double En, w, diff;
    En=0.0; w=0.0;  diff=0.0;
    for(int n=0;n<H_size;n++){
        En = Hamiltonian_BHZ_.evals_[n];
        for(int om=0;om<w_size;om++){
            w = w_min + om*dw;
            lorentzian_[n][om] = one_by_PI_ * ( (eta)/((w-En)*(w-En) + eta*eta) );
        }
    }

    //Evaluating partial FT of single particle eigenvectors, and using it to generate A(kx,ky;w) data:
    Mat_1_Complex_doub Psi_nk;
    Psi_nk.resize(H_size);
    for(int n=0; n<H_size; n++){
        Psi_nk[n] = Zero_Complex;
    }

    complex<double> Akyw_;
    double momentum_step, ky;
    int ky_ind,ky_min,ky_max;   
    if(Parameters_BHZ_.PBC_Y){
        ky_min=-ly_/2;       ky_max=ly_/2;
        momentum_step = ((2.0*PI)/(1.0*ly_));
    }
    else{
        ky_min=1;       ky_max=ly_;
        momentum_step = ((1.0*PI)/(1.0*ly_+1.0));
    }

    string file_out="Akyw_simplified.txt";
    ofstream file_Akw_out(file_out.c_str());
    file_Akw_out<<"#ky      omega       Akyw"<<endl;

    for(ky_ind=ky_min;ky_ind<=ky_max;ky_ind++){
        ky=ky_ind*momentum_step;

        //For each eigenstate |n>, computing the Fourier transformed state Psi_n(kx,ky), given by:
        //Psi_n(ky) = sum_{rx,ry} e^{-i(ky*ry)}*Psi_n(rx,ry)
        //where Psi_n(rx,ry) = evecs_[n][r]

        for(int n=0; n<H_size; n++){
            complex<double> sum_nk = Zero_Complex;

            for(int ry=0; ry<ly_; ry++){
                for(int rx=0; rx<lx_; rx++){
                    for(int orb=0; orb<N_orbs_;orb++){
                        for(int spin=0; spin<N_spin_;spin++){
                            
                            int r = Hamiltonian_BHZ_.makeIndex(rx, ry, orb, spin);
                            sum_nk += ( 1.0/sqrt(1.0*ly_) ) * ( exp(-1.0*Iota_Complex* (ky*(1.0*ry))) * Hamiltonian_BHZ_.evecs_[n][r] );
                        }
                    }
                }
            }
            Psi_nk[n] = sum_nk;
        }

        //Computing the Spectral Function at each frequency 'w' and momentum 'ky', using
        //A(ky;w) = sum_n |Psi_n(k)|^2 * lorentzian_[n][om]
        for(int om=0;om<w_size;om++){
            w = w_min + om*dw;
            Akyw_ = Zero_Complex;

            for(int n=0; n<H_size; n++){
                Akyw_ += ( conj(Psi_nk[n]) * Psi_nk[n] ) * lorentzian_[n][om];
            }
            file_Akw_out<<ky<<"   "<<w<<"      "<<Akyw_.real()<<endl;
        }

        file_Akw_out<<endl;
    }


}

/*
void Observables_BHZ::calculateBmat(){
    Mat_5_Complex_doub B_mat;
    B_mat.resize(lx_);
    for(int r1x=0;r1x<lx_;r1x++){
        B_mat[r1x].resize(lx_);

        for(int r2x=0;r2x<lx_;r2x++){
            B_mat[r1x][r2x].resize(ly_);

            for(int r1y=0;r1y<ly_;r1y++){
                B_mat[r1x][r2x][r1y].resize(ly_);

                for(int r2y=0;r2y<ly_;r2y++){
                    B_mat[r1x][r2x][r1y][r2y].resize(w_size);
                }
            }
        }
    }
    const auto& evals = Hamiltonian_BHZ_.evals_;

    //Setting the number of threads for OpenMP:
    #ifdef _OPENMP
        int threads = Parameters_BHZ_.npthreads_;
        omp_set_num_threads(threads);
        cout << "OpenMP is enabled. Number of threads: " << omp_get_max_threads() << endl;
    #else
        cout << "OpenMP is not enabled." << endl;
    #endif

    double start_time=0.0,end_time=0.0;

    #ifdef _OPENMP
        #pragma omp parallel for collapse(5) default(shared)
    #endif
    
    for(int r1x=0;r1x<lx_;r1x++){
        for(int r2x=0;r2x<lx_;r2x++){
            for(int r1y=0;r1y<ly_;r1y++){
                for(int r2y=0;r2y<ly_;r2y++){
                    for(int om=0;om<w_size;om++){
                        B_mat[r1x][r2x][r1y][r2y][om]=Zero_Complex;

                        if (om == 0) {
                            #ifdef _OPENMP
                                #pragma omp critical
                                cout << "Thread " << omp_get_thread_num() << endl;
                            #endif
                            }   
                    
                    }
                }
            }
        }
    }
    
    cout<<"B_mat initialized"<<endl;

    double w;
//    w=0.0;
    int r1,r2;
    complex<double> B_val;

    #ifdef _OPENMP
        start_time = omp_get_wtime();
        #pragma omp parallel for collapse(5) private(w, B_val, r1, r2)    
    #endif 
    for(int r1x=0;r1x<lx_;r1x++){
        for(int r2x=0;r2x<lx_;r2x++){

            for(int r1y=0;r1y<ly_;r1y++){
                for(int r2y=0;r2y<ly_;r2y++){

                    for(int om=0;om<w_size;om++){
                        w=w_min+om*dw;
                        B_val=Zero_Complex;


                        for(int orb=0;orb<N_orbs_;orb++){
                            for(int spin=0;spin<N_spin_;spin++){

                                r1 = spin * MHS + (orb + 2 * r1y + 2 * ly_ * r1x);
                                r2 = spin * MHS + (orb + 2 * r2y + 2 * ly_ * r2x);

                                for(int n=0;n<H_size;n++){
                                    double diff = w - evals[n];
                                    double denom = diff * diff + eta * eta;
                                    B_val += (one_by_PI_) * calculateAmplitude(n, r1, n, r2) * (eta / denom);

                                    //B_val += (one_by_PI_)*(calculateAmplitude(n, r1, n, r2))*
                                      //  ( (eta)/((w-Hamiltonian_BHZ_.evals_[n])*(w-Hamiltonian_BHZ_.evals_[n])+(eta*eta)) );
                                }

                            }
                        }
                        if (om == 0) {
                            #ifdef _OPENMP
                                #pragma omp critical
                                cout << "Thread " << omp_get_thread_num() << " computed B_val = " << B_val << endl;
                            #endif
                        }

                        B_mat[r1x][r2x][r1y][r2y][om] = B_val;
                    }

                }
            }
        }
    }
    
    #ifdef _OPENMP
    end_time = omp_get_wtime();
    cout << "Execution time for Bmat: " << (end_time - start_time) << " seconds" << endl;
    #endif


}


void Observables_BHZ::calculateYMomentumResolvedAkxw(){
       
    int kx_ind,ky_ind;
    int kx_min,kx_max;
    int ky_min,ky_max;
    double kx,ky,kx_step,ky_step;

    if(Parameters_BHZ_.PBC_X){
        kx_min=-lx_/2;       kx_max=lx_/2;
        kx_step = ((2.0*PI)/(1.0*lx_));
    }
    else{
        kx_min=1;       kx_max=lx_;
        kx_step = ((1.0*PI)/(1.0*lx_+1.0));
    }

    if(Parameters_BHZ_.PBC_Y){
        ky_min=-ly_/2;       ky_max=ly_/2;
        ky_step = ((2.0*PI)/(1.0*ly_));
    }
    else{
        ky_min=1;       ky_max=ly_;
        ky_step = ((1.0*PI)/(1.0*ly_+1.0));
    }


    double w;
    w=0.0; 
    complex<double> Akxw_;

    string base(".txt");
    string head("Akxw_for_ky_index_");

    for(ky_ind=ky_min;ky_ind<=ky_max;ky_ind++){
        ky=ky_ind*ky_step;

        string file_Akxw(head + to_string(ky_ind) + base);
        ofstream file_Akxw_out(file_Akxw.c_str());

        for(kx_ind=kx_min;kx_ind<=kx_max;kx_ind++){
            kx=kx_ind*kx_step;
            for(int om=0;om<w_size;om++){
                w=w_min+om*dw;
                Akxw_ = Zero_Complex;

                for(int r1x=0;r1x<lx_;r1x++){
                    for(int r2x=0;r2x<lx_;r2x++){

                        for(int r1y=0;r1y<ly_;r1y++){
                            for(int r2y=0;r2y<ly_;r2y++){

                                if(Parameters_BHZ_.PBC_X && Parameters_BHZ_.PBC_Y){
                                    Akxw_ += (1.0/(1.0*lx_*ly_))*( exp(Iota_Complex*kx*(1.0*(r1x-r2x)))*exp(Iota_Complex*ky*(1.0*(r1y-r2y)))
                                                *B_mat[r1x][r2x][r1y][r2y][om] );
                                }
                                else if(!Parameters_BHZ_.PBC_X && Parameters_BHZ_.PBC_Y){
                                    Akxw_ += (2.0/((lx_+1)*1.0*ly_))*( sin(kx*(1.0*r1x+1.0))*sin(kx*(1.0*r2x+1.0))*exp(Iota_Complex*ky*(1.0*(r1y-r2y)))
                                                *B_mat[r1x][r2x][r1y][r2y][om] );
                                }
                                else if(Parameters_BHZ_.PBC_X && !Parameters_BHZ_.PBC_Y){
                                    Akxw_ += (2.0/((ly_+1)*1.0*lx_))*( sin(ky*(1.0*r1y+1.0))*sin(ky*(1.0*r2y+1.0))*exp(Iota_Complex*kx*(1.0*(r1x-r2x)))
                                                *B_mat[r1x][r2x][r1y][r2y][om] );
                                }
                                else{
                                    Akxw_ += (4.0/((lx_+1)*(ly_+1)*1.0))*( sin(kx*(1.0*r1x+1.0))*sin(kx*(1.0*r2x+1.0))
                                                *sin(ky*(1.0*r1y+1.0))*sin(ky*(1.0*r2y+1.0))*B_mat[r1x][r2x][r1y][r2y][om] );
                                }

                            }
                        }
                    }
                }
                file_Akxw_out<<kx<<"    "<<w<<"     "<<Akxw_.real()<<endl;
            }
            file_Akxw_out<<endl;
        }
        file_Akxw_out.close();
    }

}

*/

void Observables_BHZ::calculateWaveFunctions(){

    double E_min=-0.2, E_max=0.2, eval;
    eval = 0.0;
    int state_, rs_up, rs_dn, rp_up, rp_dn;
    complex<double>val_s_up,val_s_dn, val_p_up,val_p_dn;
    
    string head = "wave_function_for_state_";
    string base = ".txt";

    for(int np=0;np<Hamiltonian_BHZ_.evals_.size();np++){
        eval = Hamiltonian_BHZ_.evals_[np];

        if(eval > E_min && eval < E_max){
            state_=np;
            val_s_up = Zero_Complex;    val_s_dn = Zero_Complex;
            val_p_up = Zero_Complex;    val_p_dn = Zero_Complex;

            string file_wave_fn_(head + to_string(state_) + base);
            ofstream file_out(file_wave_fn_.c_str());
            if (!file_out.is_open()) {
                cerr << "Error: Could not open file " << file_wave_fn_ << endl;
                continue;  // Skip to the next state if file opening fails
            }

            file_out << "#r=ry+ly*rx    phi(r,0,0)    phi(r,0,1)  phi(r,1,0)  phi(r,1,1)" << endl;

            for(int rx=0;rx<lx_;rx++){
                for(int ry=0;ry<ly_;ry++){
                    rs_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 0);
                    rs_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 0, 1);
                    rp_up = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 0);
                    rp_dn = Hamiltonian_BHZ_.makeIndex(rx, ry, 1, 1);
                                
                    val_s_up = calculateAmplitude(state_, rs_up, state_, rs_up);
                    val_s_dn = calculateAmplitude(state_, rs_dn, state_, rs_dn);
                    val_p_up = calculateAmplitude(state_, rp_up, state_, rp_up);
                    val_p_dn = calculateAmplitude(state_, rp_dn, state_, rp_dn);
                    
                    file_out<<ry+ly_*rx<<"  "<<val_s_up.real()<<"   "<<val_s_dn.real()<<"   "<<val_p_up.real()<<
                    "   "<<val_p_dn.real()<<endl;
                     
                }
            }
            file_out.close();
        }
    }

}


void Observables_BHZ::calculateEdgeConductance(){
    //G_yy = (2*PI/Ly) sum_{m,n} (f_m - f_n)/(z*z + (e_n-e_m)*(e_n-e_m))*[|<m|Jy(rx=0)|n>|^2]
    Mat_2_Complex_doub CC_left;
    CC_left.resize(H_size);
    for(int m=0;m<H_size;m++){
        CC_left[m].resize(H_size);
        for(int n=0;n<H_size;n++){
            CC_left[m][n] = Zero_Complex;
        }
    }

    double eta = 0.002;
    double mu_min = -0.2;
    double mu_max = 0.2;
    double d_mu = 0.001;
    int mu_size = (int)((mu_max-mu_min)/d_mu) + 1;
    double temp_mu;
    double e_m,e_n;
    cout<<"total-mu points="<<mu_size<<endl;

    complex<double> temp_J_aa,temp_J_ab;
    for(int m=0;m<H_size;m++){
        for(int n=0;n<H_size;n++){
            temp_J_aa = Zero_Complex;   temp_J_ab = Zero_Complex;

            for(int ry=0;ry<ly_;ry++){
                for(int orb=0;orb<N_orbs_;orb++){
                    for(int spin=0;spin<N_spin_;spin++){
                        int r = Hamiltonian_BHZ_.makeIndex(0,ry,orb,spin);
                        int r1 = Hamiltonian_BHZ_.makeIndex(0,(ry+1)%ly_,orb,spin);
                        int r2 = Hamiltonian_BHZ_.makeIndex(0,(ry+1)%ly_,1-orb,spin);
                        
                        temp_J_aa += Iota_Complex*(Hamiltonian_BHZ_.C_mat[r1][r] * calculateAmplitude(m, r1, n, r) 
                                                        - Hamiltonian_BHZ_.C_mat[r][r1] * calculateAmplitude(m, r, n, r1));
                        
                        temp_J_ab += Iota_Complex*(Hamiltonian_BHZ_.C_mat[r2][r] * calculateAmplitude(m, r2, n, r) 
                                                        - Hamiltonian_BHZ_.C_mat[r][r2] * calculateAmplitude(m, r, n, r2));
                    }
                }
            }
            CC_left[m][n] = temp_J_aa + temp_J_ab;
        }
    }

    string file_EC_yy="edge_conductance_vs_mu.txt";
    ofstream file_EC_out(file_EC_yy.c_str());

    double EC_yy;
    for(int mu_ind=0;mu_ind<mu_size;mu_ind++){
        EC_yy = 0.0;
        temp_mu = mu_min + mu_ind*d_mu;

        for(int m=0;m<H_size;m++){
            e_m = Hamiltonian_BHZ_.evals_[m];

            for(int n=0;n<H_size;n++){
                e_n = Hamiltonian_BHZ_.evals_[n];

                if(m!=n){
                    EC_yy += ((2.0*PI)/(1.0*ly_)) * (fermifunction(e_m,temp_mu)-fermifunction(e_n,temp_mu)) *
                                ( 1.0 / (eta*eta + (e_m - e_n)*(e_m - e_n)) ) * (norm(CC_left[m][n]));
                }
                
            }
        }
        file_EC_out<<temp_mu<<"   "<<EC_yy<<endl;
    }

}

#endif
