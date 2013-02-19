#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#include "healpix_base.h"
#include "add_particles.h"


//void ang2vec(double theta, double phi, vec3 *vec);

int add_particles(int Np, double *theta, double *phi, double *flux, double *radius, int Nside, double **map) {

  int i,j;

  Healpix_Base base (Nside, RING, SET_NSIDE);

  int count=0;
  int Npix_map = 12 * Nside*Nside;
  double dOmega = 4.0 * 3.1415927 / Npix_map;
  double theta0 = acos( 1.0 - dOmega/(2.0*3.141592) );

  double *weight;
  weight = (double *) malloc(Npix_map * sizeof(double));
  
  for(i=0;i<Np;i++) {

    pointing p (theta[i],phi[i]);

    vec3 vec;
    ang2vec(theta[i],phi[i],&vec);

    int pix = base.ang2pix(p);
       
    // add particle to this pixel, if it doesn't cover more than one
    if( 2.0*radius[i] < theta0 ) {
      count++;
      (*map)[pix] += flux[i];
      continue;
    }

    vector<int> pix_list;
    base.query_disc(p, 2.0*radius[i], pix_list);

    int Npix = pix_list.size();

    // if either zero or one pixel are covered by particle (this should be avoided by the above condition...)
    if(Npix < 2) {
      (*map)[pix] += flux[i];
      continue;
    }

    // get here only if the particle covers more than one pixel

    double weight_norm = 0.0;

    for(j=0;j<Npix;j++) {
      
      int this_pix = pix_list[j];
      vec3 this_vec = base.pix2vec(this_pix);
      
      double d2 = acos( dotprod(this_vec,vec) ) / radius[i];
      d2 = d2*d2;
      
      // from a Gaussian fit to the projected SPH kernel...
      weight[j] = exp(-0.5 * d2 / 0.333);      
      weight_norm += weight[j];
    
    }

    // apply weighted flux to map
    for(j=0;j<Npix;j++) {

      // first normalize weight array
      weight[j] /= weight_norm;

      (*map)[pix_list[j]] += weight[j] * flux[i];

    }  // loop over pixels covered by this particle

  }  // loop over particles

  free(weight);

  return 0;
}


void ang2vec(double theta, double phi, vec3 *vec) {

  double sz;
  double PI=M_PI;

  if( theta<0. || theta>PI) {
    fprintf(stderr, "%s (%d): theta out of range: %f\n", __FILE__, __LINE__, theta);
    exit(0);
  }

  sz = sin(theta);

  vec->x = sz * cos(phi);
  vec->y = sz * sin(phi);
  vec->z = cos(theta);

}

