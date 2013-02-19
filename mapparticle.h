#ifndef __MAP_PARTICLE_STRUCTURE__
#define __MAP_PARTICLE_STRUCTURE__
#define POS_FLOAT

#ifdef POS_FLOAT
  #define Real float
#else
  #define Real double
#endif


struct MapParticle{
	Real mass;
	Real density;
	Real hsmooth;
	Real xpos;
	Real ypos;
	Real zpos;
};



#endif

