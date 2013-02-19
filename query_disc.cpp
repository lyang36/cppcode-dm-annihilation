#include <vector>
#include "lsconstants.h"

using namespace std;


template<typename I, typename F> inline I ifloor (F arg){
	using namespace std;
	return I(floor(arg));
}

int ring_above (long nside_, double z){
	double az=abs(z);
	if (az>twothird) // polar caps
	{
		int iring = int(nside_*sqrt(3*(1-az)));
		return (z>0) ? iring : 4*nside_-iring-1;
	}
	else // ----- equatorial region ---------
	return int(nside_*(2-1.5*z));
}

void in_ring(long nside_, int iz, double phi0, double dphi,
  vector<int> &listir){
	long npface_ = nside_ * nside_;
	long ncap_   = (npface_-nside_)<<1;
	long npix_   = 12 * npface_;
	double fact2_  = 4. / (double) npix_;
	double fact1_  = (nside_<<1) * (double) fact2_;

	int nr, ir, ipix1;
	double shift=0.5;

	if (iz<nside_) // north pole
	{
		ir = iz;
		nr = ir*4;
		ipix1 = 2*ir*(ir-1);        //    lowest pixel number in the ring
	}
	else if (iz>(3*nside_)) // south pole
	{
		ir = 4*nside_ - iz;
		nr = ir*4;
		ipix1 = npix_ - 2*ir*(ir+1); // lowest pixel number in the ring
	}
	else // equatorial region
	{
		ir = iz - nside_ + 1;           //    within {1, 2*nside + 1}
		nr = nside_*4;
		if ((ir&1)==0) shift = 0;
		ipix1 = ncap_ + (ir-1)*nr; // lowest pixel number in the ring
	}

	int ipix2 = ipix1 + nr - 1;       //    highest pixel number in the ring

	// ----------- constructs the pixel list --------------
	if (dphi > (pi-1e-7))
		for (int i=ipix1; i<=ipix2; ++i) listir.push_back(i);
	else
	{
		int ip_lo = ifloor<int>(nr*inv_twopi*(phi0-dphi) - shift)+1;
		int ip_hi = ifloor<int>(nr*inv_twopi*(phi0+dphi) - shift);
		int pixnum = ip_lo+ipix1;
		if (pixnum<ipix1) pixnum += nr;
		for (int i=ip_lo; i<=ip_hi; ++i, ++pixnum)
		{
			if (pixnum>ipix2) pixnum -= nr;
				listir.push_back(pixnum);
		}
	}
}

void query_disc (long nside_, const double theta, const double phi, 
	const double radius, std::vector<int> &listpix){
	listpix.clear();
	long npface_ = nside_ * nside_;
	long ncap_   = (npface_-nside_)<<1;
	long npix_   = 12 * npface_;
	double fact2_  = 4. / (double) npix_;
	double fact1_  = (nside_<<1) * (double) fact2_;

	double dth1 = fact2_;
	double dth2 = fact1_;
	double cosang = cos(radius);

	double z0 = cos(theta);
	double xa = 1./sqrt((1-z0)*(1+z0));

	double rlat1  = theta - radius;
	double zmax = cos(rlat1);
	int irmin = ring_above (nside_, zmax)+1;

	if (rlat1<=0) // north pole in the disc
	for (int m=1; m<irmin; ++m) // rings completely in the disc
		in_ring (nside_, m, 0, pi, listpix);

	double rlat2  = theta + radius;
	double zmin = cos(rlat2);
	int irmax = ring_above (nside_, zmin);

	// ------------- loop on ring number ---------------------
	for (int iz=irmin; iz<=irmax; ++iz) // rings partially in the disc
	{
		double z;
		if (iz<nside_) // north polar cap
		z = 1.0 - iz*iz*dth1;
		else if (iz <= (3*nside_)) // tropical band + equat.
		z = (2*nside_-iz) * dth2;
		else
		z = -1.0 + (4*nside_-iz)*(4*nside_-iz)*dth1;

		// --------- phi range in the disc for each z ---------
		double x = (cosang-z*z0)*xa;
		double ysq = 1-z*z-x*x;
		//planck_assert(ysq>=0, "error in query_disc()");
		double dphi=atan2(sqrt(ysq),x);
		in_ring (nside_, iz, phi, dphi, listpix);
	}

	if (rlat2>=pi) // south pole in the disc
		for (int m=irmax+1; m<(4*nside_); ++m)  // rings completely in the disc
			in_ring (nside_, m, 0, pi, listpix);
}