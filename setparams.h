#ifndef __SETPARAMS__
#define __SETPARAMS__
using namespace std;

class SetParams{
	Info info;
	Units units;
	Natconst natconst;
	Halo halo;
	Codeunits codeunits;
	Params params;
	double rotmatrix[3][3];
	Master master;
	double randomu();
	double Dsun_in_kpc;

	/*---------------------*/
	public:
	//keywords
	double *centerpos;
	double *centvel;
	double *observerpos;
	double *particle_masses;
	double *particle_numbers;
	double *redshift;
	double *scaled;
	double *length_scale_factor;
	double *Mvir_in_Msun;
	double *M200_in_Msun;
	double *M200crit_in_Msun;
	double *Rvir_in_Mpc;
	double *R200_in_Mpc;
	double *R200crit_in_Mpc;
	double *Lbox_in_Mpc;
	double *rconverg_in_kpc;
	double *shape_radii;
	double *halo_shape;
	double *detector;
	Dm *dm;
	Map *map;
	double *analytical;
	double *other;
	double *files;
	////////////////////////////////////////


	//functions
	Master set_master(string * info_filename);
	
	void set_observer_position(Master * master, double observerpos[3], bool randompos);
	void rotation_matrix(Master * master, double align_vector[3]);
};

#endif
