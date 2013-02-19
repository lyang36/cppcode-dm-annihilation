//#include <stdio>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "cblas.h"
#include <rpc/xdr.h>
#include <rpc/types.h>
#include "tipsydefs.h"
#include "info.h"
#include "structures.h"
#include "setparams.h"
#include "readfiles.h"
#include "prepare.h"
#include "skymap.h"
#include "healpix_base.h"
#include "add_particles.h"
#include "VL2_debug.h"
#include "mapparticle.h"

//using namespace std;

int main(int argc, const char **argv){
	bool only_readwrite = false;
	int m=1;
	string basedir = "/home/lyang/files/r200";
	string basename =  "vl2b.00400.r200";
	string info_file =  "VL2_info.txt";
	string fits_file = "skymap.fits";
	string rwfile = "data.bin";
	while (m<argc)
	{
		string arg = argv[m];
		if (arg == "-base") { basedir = argv[m+1]; m+=2;}
		if (arg == "-bname") { basename = argv[m+1]; m+=2;}
		if (arg == "-info") { info_file = argv[m+1]; m+=2;}
		if (arg == "-fits") { fits_file = argv[m+1]; m+=2;}
		if (arg == "-readwrite") { only_readwrite = true; rwfile = argv[m+1]; m+=2;}
		if (arg == "-help") {
			cout << "Usage:" << endl;
			cout << "-base name: setup basedir of the data. Default: /home/lyang/files/r200" << endl;
			cout << "-bname name: setup basename of the data. Default: vl2b.00400.r200" << endl;
			cout << "-info name: setup info file path. Default: VL2_info.txt" << endl;
			cout << "-fits name: setup output fitsname. Default: skymap.fits" << endl;
			cout << "-readwrite name: read-write mode. Convert the data to map native bin data to the filename" << endl;
			cout << "-help: watching help" << endl;
			exit(0);
		}
	}
	cout << "BASE DIR: " << basedir << endl;
	cout << "BASE NAME: " << basename << endl;
	cout << "INFO FILE: " << info_file << endl;
	cout << "OUTPUT FITS: " << fits_file << endl;

/*----------------------------------------------------*/
	Prepare * prepare = new Prepare();
	Skymap * skymap = new Skymap();
	prepare->basedir = basedir;
	prepare->basename = basename;
	prepare->info_file = info_file;
	skymap->fits_filename = &fits_file;
	Master master;
	Allparts allp;
	//cout << "Good" << endl;
//	prepare->master_only = true;
	prepare->doit_VL2( master, & allp);
/*-----------------------------------------------------*/
/*
write the file tobinary data for c# reading 
*/
/*	cout << "writing data to files ... " << endl;
	int nums = prepare->Nparticle;
	cout << "There are " << nums << " particles" << endl;
	
	ofstream myFile ("data.bin", ios::out | ios::binary);;
	if(myFile.good()){
		cout << "writing mass..." << endl;
   		myFile.write ((char*)&nums, sizeof (nums));
		myFile.write ((char*)allp.mass, sizeof (double) * nums);
		cout << "writing density..." << endl;
        	myFile.write ((char*)&nums, sizeof (nums));
        	myFile.write ((char*)allp.density, sizeof (float) * nums);
        	cout << "writing hsmooth..." << endl;
        	myFile.write ((char*)&nums, sizeof (nums));
        	myFile.write ((char*)allp.hsmooth, sizeof (float) * nums);
        	cout << "writing xpos..." << endl;
        	myFile.write ((char*)&nums, sizeof (nums));
        	myFile.write ((char*)allp.xpos, sizeof (double) * nums);
        	cout << "writing ypos..." << endl;
        	myFile.write ((char*)&nums, sizeof (nums));
        	myFile.write ((char*)allp.ypos, sizeof (double) * nums);
        	cout << "writing zpos..." << endl;
        	myFile.write ((char*)&nums, sizeof (nums));
        	myFile.write ((char*)allp.zpos, sizeof (double) * nums);
		cout << "finishing ... " << endl;
		myFile.close();
	}
	exit(0);	
*/

	if(only_readwrite){
      		cout << "writing data to files ... " << endl;
       	 	int nums = prepare->Nparticle;
        	cout << "There are " << nums << " particles" << endl;
        	ofstream myFile (rwfile.c_str(), ios::out | ios::binary);;
        	MapParticle mp;
		myFile.write ((char*)&nums, sizeof (nums)); 
		if(myFile.good()){
			cout << "starting ... "<< endl;
			for(int i = 0; i < nums; i++){
				mp.mass = allp.mass[i];
				mp.density = allp.density[i];
				mp.hsmooth = allp.hsmooth[i];
				mp.xpos = allp.xpos[i];
				mp.ypos = allp.ypos[i];
				mp.zpos = allp.zpos[i];
				myFile.write ((char*)&mp, sizeof (mp));
				if( ((long)i % ((long) nums * 5l / 100l))  == 0)
					cout << " " << (100l * (long) i / (long)nums) << "% ";
					
      			}
			cout << endl;
            cout << "ending..." << endl;
            myFile.close();
		}
       	 	exit(0);
	}




//*--------------generating fits map------------------*/
	skymap->fits_filename = & fits_file;
	skymap->Np = prepare->Nparticle;
	skymap->master = &master;
	skymap->particles = &allp;
	skymap->map = &(master.map);
	skymap->creat_map();
/*----------------------------------------------------*/
	return 0;	
}
