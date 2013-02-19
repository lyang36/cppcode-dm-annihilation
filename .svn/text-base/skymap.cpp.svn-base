/**************************************************************/
/*                    Class: skymap                           */
/*	              Generate the skymap.                        */ 
/*           Author: Lin Yang      03/07/2012                 */
/**************************************************************/

#include "cblas.h"
#include <cmath>
#include <fstream>
#include <string>
#include <cstring>
#include "cblas.h"
#include <rpc/xdr.h>
#include <rpc/types.h>
#include "tipsydefs.h"
#include "info.h"
#include "structures.h"
#include "healpix_base.h"
#include "healpix_map.h"
#include "add_particles.h"
#include "arr.h"
#include "fitshandle.h"
#include "healpix_map_fitsio.h"
#include "skymap.h"
#include "VL2_debug.h"

using namespace std;;

Skymap::Skymap(){
	if(_reload == NULL)	reload = false;
	else reload = *_reload;

	if(_rotate == NULL) rotate = true;
	else rotate = *_rotate;
}

bool Skymap::creat_map(){
	//allocate memory for map
	const int LP = 10000;

	int i, j;
	//get rotation matrix
	if ( rotate ){
		cblas_dcopy(3, master->rotmatrix[0], 1, rotmatrix[0], 1);
		cblas_dcopy(3, master->rotmatrix[1], 1, rotmatrix[1], 1);
		cblas_dcopy(3, master->rotmatrix[2], 1, rotmatrix[2], 1);
	}else{
		double z[] = {0, 0, 0};
		cblas_dcopy(3, z, 1, rotmatrix[0], 1);
		cblas_dcopy(3, z, 1, rotmatrix[1], 1);
		cblas_dcopy(3, z, 1, rotmatrix[2], 1);
		rotmatrix[0][0] = 1.0;
		rotmatrix[1][1] = 1.0;
		rotmatrix[2][2] = 1.0;
	}
#ifdef _DEBUG
	//cout << "good1" <<endl;
#endif

	double * xpos = particles->xpos;
	double * ypos = particles->ypos;
	double * zpos = particles->zpos;
	double * mass = particles->mass;
	float * density = particles->density;
	float * hsmooth = particles->hsmooth;

	double * opos = master->params.opos;

/*------------------------------------------------------*/
	double distances;
/*-----------------------------------------------------*/


	double fluxes;//  = master.codeunits.annihilation_flux_to_cgs * density * mass / (4.0 * !pi * distances^2)

	double angular_radius;//  = hsmooth/distances


	long Nside = master->map.Nside;
	long Npix_in_map = master->map.Npix;
	double costheta;
	double theta;
	double phi;
	Healpix_Base base (Nside, RING, SET_NSIDE);
	double *weight;
	int count=0;
	int Npix_map = 12 * Nside*Nside;
	double dOmega = 4.0 * 3.1415927 / Npix_map;
	double theta0 = acos( 1.0 - dOmega/(2.0*3.141592) );
	double * allskymap = (double *) calloc(Npix_map, sizeof(double));
	weight = (double *) calloc(Npix_map, sizeof(double));

#ifdef _DEBUG
	cout << "good2 : " << Np << "  " << Npix_in_map << endl;
	print_out_master(master);
#endif
	int rec = Np / 50;
	int ll = 0;
	cout << "Creating map!!!" << endl;
	cout << "---10---20---30---40---50---60---70---80---90--100%\n";
	for( i = 0; i< Np; i++){

		if(i % rec == 0){
		cout << "#";
		cout.flush();
		}
		if( density[i] < 0.0){
			continue;
		}
		distances = sqrt( (xpos[i]-opos[0]) * (xpos[i]-opos[0]) + 
			(ypos[i]-opos[1]) * (ypos[i]-opos[1]) +
			(zpos[i]-opos[2]) *(zpos[i]-opos[2]) );
		fluxes = master->codeunits.annihilation_flux_to_cgs * density[i] * mass[i] / (4.0 * PI * distances * distances);

		calc_angles( xpos[i]-opos[0], ypos[i]-opos[1], zpos[i]-opos[2], distances, 
		opos, *rotmatrix, costheta, phi);
		theta = acos(costheta);
		angular_radius = hsmooth[i] / distances;
/*------------------------------------------add particles------------------------------------------------------*/
		{
			pointing p (theta,phi);

			vec3 vec;
			ang2vec(theta,phi,&vec);

			int pix = base.ang2pix(p);
			if( 2.0*angular_radius < theta0 ) {
				count++;
				allskymap[pix] += fluxes;
				continue;
			}

			vector<int> pix_list;
			base.query_disc(p, 2.0*angular_radius, pix_list);

			int Npix = pix_list.size();

			// if either zero or one pixel are covered by particle (this should be avoided by the above condition...)
			if(Npix < 2) {
				allskymap[pix] += fluxes;
				continue;
			}

			// get here only if the particle covers more than one pixel
			double weight_norm = 0.0;

			for(j=0;j<Npix;j++) {
				int this_pix = pix_list[j];
				vec3 this_vec = base.pix2vec(this_pix);

				double d2 = acos( dotprod(this_vec,vec) ) / angular_radius;
				d2 = d2*d2;
				weight[j] = exp(-0.5 * d2 / 0.333);      
				weight_norm += weight[j];		    
			}

			// apply weighted flux to map
			for(j=0;j<Npix;j++) {
				// first normalize weight array
				weight[j] /= weight_norm;
				allskymap[pix_list[j]] += weight[j] * fluxes;

			}  // loop over pixels covered by this particle


		}
/*****************************************************************************************************************/
	}

	cout << endl;
#ifdef _DEBUG
	//cout << "good3" <<endl;
	print_out_master(master);
#endif

	//divide by solid angle of map's pixels
	//conversion of allsky from g^2 cm^-5 to Gev^2 cm^-6 kpc

	double unit_factor = pow(pow((master -> natconst.c_in_cgs), 2) /
	(master->units.eV_in_cgs * 1.0e9), 2) / (master->units.Mpc_in_cgs * 1.0e-3);
	cblas_dscal( Npix_in_map, unit_factor / master->map.dOmega, allskymap, 1);

#ifdef _DEBUG
	cout << "flux" << fluxes << endl;
	cout << "unit_factor: " << unit_factor << endl;
	cout << "dOmega: " << master->map.dOmega << endl;
	cout << "All Skymap:" <<endl;
	cout << "skymap[0]: " << allskymap[0] << endl;
	cout << "skymap[100]: " << allskymap[100] << endl;
	cout << "skymap[10000]: " << allskymap[10000] << endl;
#endif
	arr<double> maparr (allskymap, Npix_in_map);
	Healpix_Map<double> outputmap (maparr, NEST);
	if( fits_filename != NULL ){
		ifstream ifile(fits_filename->c_str());
		if (ifile) {
			// The file exists, and is open for input
			cout << "File exists! Owerite?" << endl;
			char t;
			cin >> t;
			if( t == 'y'){
				remove( fits_filename->c_str() );
				fitshandle fitswriter;
				fitswriter.create(*fits_filename);
				write_Healpix_map_to_fits(fitswriter, outputmap, TDOUBLE );
				fitswriter.close();
			}
		}else{
			cout << "Creating fits!..." << endl;
			fitshandle fitswriter;
			fitswriter.create(*fits_filename);
			write_Healpix_map_to_fits(fitswriter, outputmap, TDOUBLE );
			fitswriter.close();
		}
	}
	free(allskymap);
	free(weight);

}


void Skymap::calc_angles( double xpos, double ypos, double zpos, double &distances, 
						 double * opos, double *rotmatrix, double & costheta, double &phi){
	double vec[] = { xpos, ypos, zpos};
	double temp[] = {0, 0, 0};
	cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, rotmatrix, 3, vec, 1, 1.0, temp, 1 );
	
	double xcoord = vec[0];
	double ycoord = vec[1];
	double zcoord = vec[2];
  
	costheta = zcoord / distances;

	phi = atan2( ycoord, xcoord );
	
	if( phi < 0 ){
		phi += 2.0 * PI;
	}

	//a few adjustments...
  
	//one more rotation by 180 degrees...
	phi -= PI;
  
	//force phi to lie between 0 and 2*pi  
	if( phi < 0 ){
		phi = 2.0 * PI + phi;
	}
  
	if( phi > 2 * PI){
		phi = phi - 2.0 * PI;
	}
}

