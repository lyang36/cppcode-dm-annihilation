#ifndef __INFO__
#define __INFO__
// Info from info file
typedef struct{
	double Omega_M;
	double h100;
	double n_s;
	double sigma_8;
	double Lbox_in_Mpc;
	double particle_masses[10];
	long particle_numbers[10];
	double centerpos[3];
	double centervel[3];
	double Mvir_in_Msun;
	double M200_in_Msun;
	double M200crit_in_Msun;
	double Rvir_in_Mpc;
	double R200_in_Mpc;
	double R200crit_in_Mpc;
	double Vmax_in_kms;
	double RVmax_in_kpc;
	double rconverg_in_kpc;
	double params_NFW[3];
	double params_GNFW[3];
	double params_Einasto[3];
}Info;
#endif
