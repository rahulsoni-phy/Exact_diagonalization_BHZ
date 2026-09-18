#ifndef Parameters_BHZ_class
#define Parameters_BHZ_class
#include <iostream>
#include <fstream>
#include "tensors.h"

/*
//this function will convert a string to a boolean in its different forms:
bool string_to_bool(string& str) {
	return (str == "True" || str == "true" || str == "ON" || str == "on");
}

//this function will remove leading and trailing whitespace characters like " " or "\t"
std::string clip(std::string& str) {
    size_t start = str.find_first_not_of(" \t");
    size_t end = str.find_last_not_of(" \t");
    return (start == std::string::npos) ? "" : str.substr(start, end - start + 1);
}
*/
//Above functions will be implemented in future!!

//this class will read from the **input file**
class Parameters_BHZ{

public:
	// Declaring system params:---->
	bool PBC_X, PBC_Y;
	bool TBC_Y;
	double Phi_y;
	int Lx_, Ly_, N_Orbs, W_;
	int Total_Sites, Total_Cells, Ham_Size;

	// Declaring model params:----->
	double A_val, B_val, M_val, D_val, C_val;
	double Vo_val,fill_;
	double Hfield, Pinfield;
	int total_particles_;
	double mu_fixed, temp_, beta;
	int numberofSoftedges_;
	bool canonical_;

	double helical_edge_pot;

	// Observables details:--->
	bool get_dos, get_ldoe;
	bool get_sc, get_Akxw, get_ldos;
	bool get_Akyw, get_ky_resolved_Akxw;
	bool save_evecs;
	double w_min, w_max, dw_, eta_;

	//----------------------->
	bool get_wave_fn;
	int state_, npthreads_;

	void Initialize(std::string input_file);
	double matchstring(std::string file, std::string match);
	std::string matchstring2(std::string file, std::string match);
};

void Parameters_BHZ::Initialize(string input_file){

	cout << "-------------------------------------------\n"
		 << "Reading the input file = " << input_file << "\n"
		 << "-------------------------------------------\n"
		 << endl;
	//------------------------------------------------------------------------//
	// Reading and setting boundary conditions
	string PBC_X_string, PBC_Y_string, TBC_Y_string;
	string pbc_x_out, pbc_y_out;

	PBC_X_string = matchstring2(input_file, "PBC_X");
	if (PBC_X_string == "True")
	{	PBC_X = true;		pbc_x_out = "PBC";		}
	else if(PBC_X_string == "False")
	{	PBC_X = false;		pbc_x_out = "OBC";		}

	TBC_Y_string = matchstring2(input_file, "TBC_Y");

	PBC_Y_string = matchstring2(input_file, "PBC_Y");
	if (PBC_Y_string == "True")
	{	PBC_Y = true;		pbc_y_out = "PBC";		
		if(TBC_Y_string == "True"){
			TBC_Y = true;
			Phi_y = matchstring(input_file, "Flux_Y");
			cout<<"Flux_Y="<<Phi_y<<endl;
			save_evecs = true;
		}
		else{
			TBC_Y = false;
			save_evecs = false;
		}
	}
	else if(PBC_Y_string=="False")
	{	PBC_Y = false;		pbc_y_out = "OBC";		}


	numberofSoftedges_ = int(matchstring(input_file, "Number_of_Softwalls"));

	Lx_ = int(matchstring(input_file, "Cells_X"));
	Ly_ = int(matchstring(input_file, "Cells_Y"));
	N_Orbs = int(matchstring(input_file, "Total_Orbs"));

	Total_Cells = Lx_ * Ly_;
	Total_Sites = N_Orbs * Total_Cells;
	Ham_Size = 2 * Total_Sites;

	cout << "Boundary conditions = " << pbc_x_out << "x" << pbc_y_out << "\n"
		 << "Total size of the Hamiltonian = " << Ham_Size << "x" << Ham_Size << endl;
	//-----------------------------------------------------------------------//

	A_val = matchstring(input_file, "Inter_Orb_Hopping_A");
	B_val = matchstring(input_file, "Intra_Orb_Hopping_B");
	D_val = matchstring(input_file, "Intra_Orb_Hopping_D");
	M_val = matchstring(input_file, "Gap_Parameter_M");
	C_val = matchstring(input_file, "Gap_Parameter_C");
	Hfield = matchstring(input_file, "Mag_Field");
	Pinfield = matchstring(input_file, "Pinning");

	helical_edge_pot = matchstring(input_file, "Helical_Edge_Pot");

	string Ensemble_str, Ens_out;
	Ensemble_str = matchstring2(input_file, "Ensemble");
	if (Ensemble_str == "CE")
	{	canonical_ = true;		Ens_out = "Canonical";		}
	else
	{	canonical_ = false;		Ens_out = "Grand-Canonical";}

	fill_ = matchstring(input_file, "Filling");
	temp_ = matchstring(input_file, "Temperature");
	beta = 1.0 / (1.0 * temp_);

	W_ = int(matchstring(input_file, "Softwall_Width"));
	Vo_val = matchstring(input_file, "Confining_Potential_Vo");
	mu_fixed = matchstring(input_file, "Fixed_mu");

	cout << "(A, B, D, M) = (" << A_val << " , " << B_val << " , " << D_val << " , " << M_val << ")" << "\n"
		 << "Performing " << Ens_out << "-ensemble" << endl;
	if (canonical_ == true)
	{
		total_particles_ = (int) (fill_*Ham_Size + 0.05);
		cout << "filling = " << fill_ << endl;
	}
	if (canonical_ == false)
	{
		cout << "Confining-width = " << W_ << endl;
		cout << "Confining-potential = " << Vo_val << endl;
	}
	//----------------------------------------------------------------------//

	string ldoe_str, dos_str, sc_str, ldos_str;
	string Akxw_str, Akyw_str, ky_res_Akxw_str;
	ldoe_str = matchstring2(input_file, "Calculate_LDOE");
	if (ldoe_str == "True")
	{	
		get_ldoe = true;
		cout << "Measuring LDOE" << endl;		
	}
	else {	get_ldoe = false;	}

	dos_str = matchstring2(input_file, "Calculate_DOS");
	if (dos_str == "True")
	{
		get_dos = true;
		cout << "Measuring DOS" << endl;
	}
	else {	get_dos = false;	}

	sc_str = matchstring2(input_file, "Calculate_SC");
	if (sc_str == "True")
	{
		get_sc = true;
		cout << "Measuring Spin Currents" << endl;
	}
	else {	get_sc = false;		}

	Akxw_str = matchstring2(input_file, "Spectral_Akxw");
	if (Akxw_str == "True")
	{
		get_Akxw = true;
		cout << "Measuing spectral function A(kx,w)" << endl;
	}
	else {	get_Akxw = false;	}

	Akyw_str = matchstring2(input_file, "Spectral_Akyw");
	if (Akyw_str == "True")
	{
		get_Akyw = true;
		cout << "Measuing spectral function A(ky,w)" << endl;
	}
	else {	get_Akyw = false;	}

	ky_res_Akxw_str = matchstring2(input_file, "Spectral_ky_resolved_Akxw");
	if (ky_res_Akxw_str == "True")
	{
		get_ky_resolved_Akxw = true;
		cout << "Measuing mom.-y resolved spectral function A(kx,w)" << endl;
	}
	else {	get_ky_resolved_Akxw = false;	}

	ldos_str = matchstring2(input_file, "Calculate_LDOS_along_x");
	if (ldos_str == "True")
	{	
		get_ldos = true;
		cout << "Measuring LDOS along-x" << endl;
	}
	else {	get_ldos = false;	}

	w_min = matchstring(input_file, "omega_min");
	w_max = matchstring(input_file, "omega_max");
	dw_ = matchstring(input_file, "d_omega");
	eta_ = matchstring(input_file, "broadening");
	//----------------------------------------------------------------------//

	string wave_fn_str;
	wave_fn_str = matchstring2(input_file, "Calculate_wave_fn");
	if (wave_fn_str == "True")
	{
		get_wave_fn = true;
		state_ = int(matchstring(input_file, "State_val"));
		cout << "Measuring wave functions" << endl;
	}
	else {	get_wave_fn = false;	}

	npthreads_ = int(matchstring(input_file, "Threads"));
	if (npthreads_ > 1)
	{
		cout << "Total threads= " << npthreads_ << endl;
	}

	cout << "-------------------------------------------\n"
		 << "Finish Reading the input file" << "\n"
		 << "-------------------------------------------\n"
		 << endl;
}

double Parameters_BHZ::matchstring(string file, string match){

	std::string test, line;
	std::ifstream inputfile(file);
	double amount;

	bool pass = false;

	while (std::getline(inputfile, line)){
		std::istringstream iss(line);

		if (std::getline(iss, test, '=') && !pass){
			if (iss >> amount && test == match){
				pass = true;
			}
			if (pass){
				break;
			}
		}
	}

	if (!pass){
		throw std::invalid_argument("Missing argument in the input file: " + match);
	}

	return amount;
}

string Parameters_BHZ::matchstring2(string file, string match){
	std::string line;
    std::ifstream inputfile(file);
    std::string amount;
    int offset;

	/* Referenced from https://cplusplus.com/forum/beginner/121556/
			 * https://stackoverflow.com/questions/12463750/ */
	if(inputfile.is_open()) {
        while(std::getline(inputfile, line)) {
            if((offset = line.find(match, 0)) != std::string::npos) {
                amount = line.substr(offset + match.length() + 1);
                break; // Break early if the match is found
            }
        }
        inputfile.close();
    } else {
        std::cerr << "Unable to open input file while in the Parameter class." << std::endl;
        return "";
    }

//    std::cout << match << " = " << amount << std::endl;
    return amount;
}


#endif
