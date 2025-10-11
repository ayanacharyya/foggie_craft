import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pickle as pkl
from collections import namedtuple
from globalpars import *

#	Auxiliary function to calculate DM related quantities

#	A named tuple to store various parameters related to radial ne profile
neradial		=	namedtuple('neradial',['logsm','logsfr','redshift','theta0','phi0','radkpc','theta','phi','inclination','neincrad'])




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
    



def radialexp (x, x0, a0):
#   Return a radial exponential
#	Arguments:	Radius, scale radius, normalization

	radexp	= a0 * np.exp(-x/x0)

	return (radexp)
#	------------------------------------------------------------------------------------------------------





















