#ifndef __STRUCTURES__
#define __STRUCTURES__

#ifdef POS_FLOAT
  #define Real float
#else
  #define Real double
#endif

using namespace std;

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

const double PI  = M_PI; 

typedef struct{
	static const double Msun_in_cgs = 1.9889212e+33;
	static const double Mpc_in_cgs = 3.08568025e+24;
	static const double eV_in_cgs = 1.602e-12;
} Units;

typedef struct{
	double h100;
	static const double G_in_cgs = 6.67259e-8;
	double rho_crit_in_cgs;
	double Omega_m;
	double Delta_vir;
        double Rvir_MW_in_Mpc;// = 0.275;
	static const double c_in_cgs = 2.99792458e10;
} Natconst;

typedef struct{
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
	double shape_r;
	double shape;
} Halo;

typedef struct{
	double M_in_GeV;// = 46.0;
	double M_in_cgs;// = 0.0;
	double sigma_v_in_cgs;// = 5.0e-26;
	
} Dm;

typedef struct{
	double mass_to_cgs;
	double mass_to_Msun;
	double length_to_Mpc;
	double length_to_cgs;
	double time_to_cgs;
	double velocity_to_cgs;
	double density_to_cgs;
	double annihilation_flux_to_cgs;
} Codeunits;

typedef struct{
	double z;
	double cpos[3];
	double cvel[3];
	double opos[3];
	double otheta;
	double ophi;
	double Lbox_in_Mpc;
	long particle_numbers[10];
	double particle_masses[10];
	double particle_masses_in_Msun[10];
} Params;


typedef tipsy_header Tipsyheader;

typedef struct{
	float mass;
	float density;
	float hsmooth;
	double xpos;
	double ypos;
	double zpos;
} Particle;

typedef struct{
	string projection;
	long Nside;
	long Npix;
	double dOmega;
	double theta0;
}Map;

typedef struct{
	double * mass;
	float * density;
	float * hsmooth;
	double * xpos;
	double * ypos;
	double * zpos;
	double * distance;
} Allparts;

typedef struct{
	Natconst natconst;
	Dm dm;
	Halo halo;
	Params params;
	Units units;
	Codeunits codeunits;
	double detector;
	Map map;
	double files;
	double analytical;
	double other;
	double (*rotmatrix)[3];//[3][3];	
} Master;

#endif
