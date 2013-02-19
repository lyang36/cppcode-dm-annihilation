/*********************************************************************/
/*              class: Prepare	                                     */
/*                 Lin Yang                                          */
/*                03/07/2012                                         */
/*  Read the particles and setup the structure "master" and particles*/
/*********************************************************************/

#include <cmath>
#include <cstring>
#include <iostream>
#include "cblas.h"
#include <rpc/xdr.h>
#include <rpc/types.h>
#include "tipsydefs.h"
#include "info.h"
#include "structures.h"
#include "setparams.h"
#include "readfiles.h"
#include "prepare.h"
#include "VL2_debug.h"


using namespace std;

Prepare::Prepare(){
	/*Setting up parameters*/
	if(_Nside == NULL){
		Nside =(long)pow(2.0, 9.0);
	}else{
		Nside = *_Nside;
	}
	if(_master_only == NULL){
		master_only = false;
	}else{
		master_only = *_master_only;
	}
	if(_reload == NULL){
		reload = false;
	}else{
		reload = *_reload;
	}

}	


void Prepare::doit_VL2(Master & master, Allparts * allp){
/*	info_file = "VL2_info.txt";
	basedir = "./r200";
	basename = "vl2b.00400.r200";*/
	/*-----------------healpix-----------------------------------------------------------------*/

	Npix = 12 * (Nside * Nside);

	dOmega = 4.0 * PI / Npix;
	theta0 = acos(1.0 - dOmega / (2.0 * PI));
	
	map.projection = "mollweide";
	map.Nside = Nside;
	map.Npix = Npix;
	map.dOmega = dOmega;
	map.theta0 = theta0;
	print_out_map( &map );
	
	/*-----------------------------------------------------------------------------------------*/
	/*--------------set master structure-------------------------------------------------------*/
	setparams = new SetParams();
	readfiles = new ReadFiles();
	double observerpos_t[] = {0.36116062, 0.21047031, 0.0069339937};
	double align_vector_t[] = {-0.0265062, -0.285643, -0.957970};

	//master = new Master;
	setparams -> map = & map;
	
	double shapes = 0;
	setparams -> halo_shape = & shapes;
	master = (*setparams).set_master(& info_file);
	(*setparams).set_observer_position( & master, observerpos_t, false);
	
	(*setparams).rotation_matrix( & master, align_vector_t);
	//cout <<"MASTER:   "<< master -> params.Lbox_in_Mpc << endl;
	print_out_master( & master );
	if(master_only){
		return;
	}
	/*-----------------------------------------------------------------------------------------*/

	/*--------------------read data------------------------------------------------------------*/
	string filename;
	filename = basedir + "/" + basename;

	readfiles->filename = &filename;
	readfiles->tipsyheader = & header;
	readfiles->native = false;
	readfiles->redshift = & redshift;
	readfiles->Nparticles = & Nparticle;
	/*--------------------read header----------------------------------------------------------*/
	(*readfiles).read_tipsyheader();
	Nparticle = header.ndark;

#ifdef _DEBUG
	Nparticle = 1000000;
	header.ndark = 1000000;
#endif

	if(redshift < 1e-3){
		redshift = 0.0;
	}	  
	
	( master).params.z = redshift;
	


#ifdef _DEBUG
	cout << "---------TIPSYHEADER----------------------" << endl;
	print_out_tipsyheader(&header);
	cout << "---------NEW MASTER----------------------" << endl;
	print_out_master(& master);

#endif
	/*--------------------read particles-------------------------------------------------------*/
	cout << "reading particles ..." << endl;
	//read the mass, xpos, ypos, and zpos to allp
	readfiles->particles = allp;
	readfiles->mass = true;
	(*readfiles).read_particles();

	/*--------------------read density----------------------------------------------------------*/
	cout << "reading densities ..." << endl;
	//filename, ".den32";
	filename = basedir + "/" + basename + ".den32";
	readfiles->filename = &filename;
	(*readfiles).read_scalar(allp->density );

	/*--------------------read hsmooth data------------------------------------------------------*/
	cout << "reading hsmooth ..." << endl;
    //strcat(filename, ".hsm32");
	filename = basedir + "/" + basename + ".hsm32";
	readfiles->filename = &filename;
    (*readfiles).read_scalar(allp->hsmooth);

#ifdef _DEBUG
	cout << "DATA:" << endl;
	cout << "POS:" << endl;
	cout << allp->xpos[0] << " " <<  allp->ypos[0] << " " <<allp->zpos[0] << endl;
	cout << allp->xpos[100] << " " <<  allp->ypos[100] << " " <<allp->zpos[100] << endl;
	cout << allp->xpos[10000] << " " <<  allp->ypos[10000] << " " <<allp->zpos[10000] << endl;
	cout << "MASS, DENS, HSMOOTH" << endl;
	cout << allp->mass[0] << " "  <<allp->density[0] << " "  <<allp->hsmooth[0] << endl;
	cout << allp->mass[100] << " "  <<allp->density[100] << " "  <<allp->hsmooth[100] << endl;
	cout << allp->mass[10000] << " "  <<allp->density[10000] << " "  <<allp->hsmooth[10000] << endl;
#endif

	cout << "setting density to -1 for all lower resolution particles..." << endl;
	double hires_particle_mass = min((master).params.particle_masses, 3);
	long finds = findwhere(allp->mass, hires_particle_mass * 1.1, allp->density, -1.0, Nparticle);
#ifdef _DEBUG
	cout <<hires_particle_mass << endl;
#endif
	cout << "Found " << finds << " High resolution particles with mass " << hires_particle_mass * master.codeunits.mass_to_Msun << " M_sun (" << hires_particle_mass << " in code unites)." << endl;
}


/*--------------------these functions is like where in IDL-----------------------------------------*/
long Prepare::findwhere(double comtar[], double lesscompare, double editar[], double falsval, int size){
	long finds = 0;
	long index = 0 ;

	for(index = 0; index < size; index ++){
		if(comtar[ index ] < lesscompare){
			finds ++;
		}
		else{
			editar[ index ] = falsval;
		}	
	}
	return finds;
}

/*--------------------these functions is like where in IDL-----------------------------------------*/
long Prepare::findwhere(double comtar[], double lesscompare, float editar[], double falsval, int size){
	long finds = 0;
	long index = 0 ;

	for(index = 0; index < size; index ++){
		if(comtar[ index ] < lesscompare){
			finds ++;
		}
		else{
			editar[ index ] = falsval;
		}	
	}
	return finds;
}
/*--------------------these functions is like min in IDL-------------------------------------------*/
double Prepare::min(double ar[], int size){
	long index;
	double minimum = ar[ 0 ];
	for(index = 0; index < size; index ++){
		if (ar[ index ]< minimum)
			minimum = ar[ index ];
	} 
	return minimum;
}
