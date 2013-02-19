#include <iostream>

using namespace std;

int add_particles(int Np, double *theta, double *phi, double *flux, double *radius, int Nside, double **map);

extern "C" int idl_wrapper(int argc, void *argv[]) {
  int Np, Nside;
  double *theta, *phi, *flux, *radius;
  double *map;


  /* Insure that the correct number of arguments were passed in */
  if(argc != 7) return 0;


  Np = *((int *) argv[0]);

  theta = (double *) argv[1];
  phi = (double *) argv[2];
  flux = (double *) argv[3];
  radius = (double *) argv[4];

  Nside = *((int *) argv[5]);

  map = (double *) argv[6];

  add_particles(Np, theta, phi, flux, radius, Nside, &map);

  return 0;
}
