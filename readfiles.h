#ifndef __READFILES__
#define __READFILES__
using namespace std;

class ReadFiles{
	
	public:
		//keywords
		//read tipsyheader
		string * filename;
		Tipsyheader * tipsyheader;
		bool native;
		long * Nparticles;
		double * redshift;

		//readparticles
		Allparts *particles;
		bool mass;
		double * value;
        bool doublepos;
		bool relative_positions;
		double centerpos[3];


		//read_scarlar

		//function
		void read_tipsyheader();
		void read_particles();
		void read_scalar(float * &); 
};

#endif
