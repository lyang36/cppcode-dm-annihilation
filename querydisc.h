#ifndef __QUERY__DISC__
#define __QUERY__DISC__

 /*! Returns the numbers of all pixels that lie at least partially within
\a radius of \a (theta, phi)  in \a listpix. It may also return a few pixels
	which do not lie in the disk at all.
\param dir the angular coordinates of the disc center
\param radius the radius (in radians) of the disc
\param listpix a vector containing the numbers of all pixels within
	the disc
\note This method works in both RING and NEST schemes, but is
	considerably faster in the RING scheme. */
void query_disc (long nside_, const double theta, const double phi, 
	const double radius, std::vector<int> &listpix);

#endif