#ifndef __SKYMAP__
#define __SKYMAP__

class Skymap{
private:
	double rotmatrix[3][3];
	bool rotate;
    bool reload;
public:
	long Np;
	Master * master;
	Allparts * particles;
	string * fits_filename;
    string * png_filename;
    double * allskymap;
	Map * map;
	Skymap();
	bool creat_map();
	/*parameters*/
	bool *_rotate;
	bool *_reload;

	void calc_angles( double xpos, double ypos, double zpos, double &distances, 
						 double * opos, double *rotmatrix, double & costheta, double &phi);

};


#endif
