from craft_utils import *

#	Auxiliary function to calculate DM related quantities


#	----------------------------------------------------------------------------------------------------------
def neprofile (necube, dkpc, theta, phi, cenpx, radius):
	#	Calculate radial electron density profiles along theta, phi	
	
	nepro	=	np.zeros(len(radius), dtype=float)
	for i in range(0, len(radius)):
		nepro[i]= necube[ int(cenpx[0]+radius[i] * np.cos(theta) * np.cos(phi)), int(cenpx[1]+radius[i] * np.cos(theta) * np.sin(phi)), int(cenpx[2] + radius[i] * np.sin(theta))]
	
	return (nepro)
#	------------------------------------------------------------------------------------------------------




def inclinvec (vecomps, theta0, phi0):
#   Return inclination (in degrees) of a vector w.r.t. z axis
#	Arguments:	A 3-component vector, orientation of the galaxy

	dotpabs	= np.abs(vecomps[0] * np.sin(theta0) * np.cos(phi0) + vecomps[1] * np.sin(theta0) * np.sin(phi0) + vecomps[2] * np.cos(theta0))
	
	inclin  = np.rad2deg(np.arccos(dotpabs / np.sqrt(np.sum(vecomps**2))))

	return (inclin)
#	------------------------------------------------------------------------------------------------------




def impactfac (pt1, pt2):
#   Return the impact factor for a given LoS in units of "pixels"
#	Arguments:	Two fixed points on the LoS
    
    crosq	= (pt1[0]*pt2[1]-pt1[1]*pt2[0])**2 +(pt1[0]*pt2[2]-pt1[2]*pt2[0])**2 + (pt1[2]*pt2[1]-pt1[1]*pt2[2])**2
    modsq	= (pt1[0]-pt2[0])**2 + (pt1[1]-pt2[1])**2 + (pt1[2]-pt2[2])**2
    
    impactf = np.sqrt(crosq / modsq) 
    
    return (impactf)
#	------------------------------------------------------------------------------------------------------
    


def logradialexp (x, x0, a0):
#   Return a radial exponential
#	Arguments:	Radius, scale radius, normalization

	radexp	= np.log10(a0 * np.exp(-x/x0))

	return (radexp)
#	------------------------------------------------------------------------------------------------------



def logradialgaus (x, x0, a0):
#   Return a radial exponential
#	Arguments:	Radius, scale radius, normalization

	radexp	= np.log10(a0 * np.exp(-(x/x0)**2))

	return (radexp)
#	------------------------------------------------------------------------------------------------------




def logradialexp3 (x, x0, a0):
#   Return a radial exponential
#	Arguments:	Radius, scale radius, normalization

	radexp	= np.log10(a0 * np.exp(-(x/x0)**0.33))

	return (radexp)
#	------------------------------------------------------------------------------------------------------



def logbeselk (x, x0, a0):
#   Return a radial Bessel k
#	Arguments:	Radius, scale radius, normalization

	radexp	= np.log10(a0 * special.kn(0.5,x/x0))

	return (radexp)
#	------------------------------------------------------------------------------------------------------


def radialexp (x, x0, a0):
#   Return a radial exponential
#	Arguments:	Radius, scale radius, normalization

	radexp	= a0 * np.exp(-x/x0)

	return (radexp)
#	------------------------------------------------------------------------------------------------------




def radialexpower (x, x0, a0, power):
#   Return a radial powered exponential
#	Arguments:	Radius, normalization, exponent

	radexp	= np.log(a0) - ((x/x0)**power)

	return (radexp)
#	------------------------------------------------------------------------------------------------------





def distfrmajorax (theta0, phi0, pt1, pt2):
#   Return the projected distance of the LoS from the apparent major axis
#	Arguments:	LoS vector, Normal vector
    
	normvec		= np.array([ np.sin(theta0) * np.cos(phi0), np.sin(theta0) * np.sin(phi0), np.cos(theta0)])

	losvec		= pt2 - pt1

	majaxpt		= np.array([0.0, 0.0, 0.0])
	majaxvec	= np.cross(losvec, normvec)

	crpdct		= np.cross(majaxvec, losvec)
	mindist		= np.abs( np.dot(crpdct, (pt1 - majaxpt)) ) / np.sqrt(np.dot(crpdct, crpdct))
	
	return (mindist)
#	------------------------------------------------------------------------------------------------------


def linearxy (xy, a, b, c):
#	Retruns a linear function of two independent variables
# 	z	=	ax + by + c	

	z	= a*xy[0] + b*xy[1] + c 

	return (z)
#	------------------------------------------------------------------------------------------------------















