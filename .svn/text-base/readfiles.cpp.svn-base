/**************************************************************/
//					Class: ReadFiles
//			Read particle data from files
//
/**************************************************************/

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include "cblas.h"
#include <rpc/xdr.h>
#include <rpc/types.h>
#include "tipsydefs.h"
#include "info.h"
#include "structures.h"
#include "readfiles.h"
#include "VL2_debug.h"

using namespace std;

void ReadFiles::read_tipsyheader(){
	char fi[255];
	strcpy(fi, filename->c_str());
	read_tipsy_header(fi, (tipsy_header *) tipsyheader);
	*Nparticles = tipsyheader -> ndark;
	*redshift = 1.0 / tipsyheader -> time - 1.0;

}


void ReadFiles::read_particles(){
	long i,j;
	Pdm dp;
	FILE * fp;
	int status;
	XDR xdrs;
	int Nparticles = 0;
	double redshift = 0;	
	char fi[255];

	if(relative_positions){
		if(centerpos == NULL){
			cout << "No center pos!!" << endl;
			exit(1);
		}
		cout << "WARNING: relative velocity is not working yet!" << endl;
	}
	//keywords
	//positions keyword overrides the {xyz}pos keywords
 
	strcpy(fi, filename->c_str());
	fp=fopen(fi,"r");
	if(!fp) {
		fprintf(stderr,"Problem opening file: %s\n",filename);
		exit(1);
		return ;
	}
	
	Nparticles = tipsyheader -> ndark;
	redshift = 1.0 / (tipsyheader -> time) - 1.0;
	xdrstdio_create(&xdrs,fp,XDR_DECODE);
	//read header
	status = xdr_header(&xdrs, tipsyheader);
	//buid all particls
	particles->mass = new double[Nparticles];
	particles->xpos = new double[Nparticles];
	particles->ypos = new double[Nparticles];
	particles->zpos = new double[Nparticles];

	int rec = Nparticles/50;
	cout << "N = " << Nparticles <<  endl;
	cout << "---10---20---30---40---50---60---70---80---90--100%\n";
	for(i=0; i < Nparticles; i++) {
		status = xdr_dark(&xdrs,&dp);
		if (status != 1) {
			fprintf(stderr,"Error reading dark particle from input file.\n");
			exit(1);
		}
		particles->xpos[i] =  dp.pos[0];
		particles->ypos[i] =  dp.pos[1];
		particles->zpos[i] =  dp.pos[2];
		particles->mass[i] =  dp.mass;
		if(i % rec == 0){
			cout << "#";
			cout.flush();
		}
	}
	cout << endl;
	xdr_destroy(&xdrs);
	fclose(fp);
}


void ReadFiles::read_scalar(float * &scalar){
	scalar = new float[ *Nparticles ];
	int np;
    ifstream myFile ( filename->	c_str(), ios::in | ios::binary);
	myFile.read( reinterpret_cast<char*>( &np ), sizeof np );

#ifdef _DEBUG
	np = * Nparticles;
#endif

    if ( np != (* Nparticles)) {
		cout << "ERROR! " << np << "---" << *Nparticles << endl;
		exit(1);
    }
#ifdef _DEBUG
	cout << "NP: " <<  np << "---" << *Nparticles << endl;
#endif

    if (!myFile.read (reinterpret_cast<char*>( scalar ), (sizeof(float)) * np )) {
		cout << "READ ERROR!"  << endl;
		exit(1);
    }
#ifdef _DEBUG
	cout << "Scale: " <<  scalar[0] << endl;
	cout << "Scale: " <<  scalar[100] << endl;
	cout << "Scale: " <<  scalar[10000] << endl;
#endif
	myFile.close();

}

