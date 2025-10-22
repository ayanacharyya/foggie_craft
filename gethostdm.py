#
#	Script for estimating FRB host galaxy DM from simulated electron density cubes
#
#								AB, August 2024

#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
import pickle as pkl
from collections import namedtuple
from globalpars import *
from nefns import *
from plotdm import *

def print_instructions():

	#	Print instructions to terminal
	
	print("\n            You probably need some assistance here!\n")
	print("\n Arguments are       --- <FITS file name> <mode> <nfixpts> <scale/kpc>\n")
	print(" Supported Modes are --- profile        (calculate electron density profiles)")
	print("                     --- losdm          (calculate LoS DMs)")
	print("                     --- pltdm          (Plot LoS DMs)")
	print("                     --- dmscat         (Plot LoS DMs)")
	
	print("\n            Now let's try again!\n")
	
	return(0)

#	--------------------------	Read inputs	-------------------------------
if(len(sys.argv)<3):
	print_instructions()
	sys.exit()

fitsname	=	sys.argv[1]					#	FITS file name
exmode		=	sys.argv[2]					#	What to do	
nfixpts     =   int(sys.argv[3])            #   Number of fixed points on each face to simulate LoSs
scalekpc	=	float(sys.argv[4])			#	Scale radius in kpc

incranges	=	np.array([[0,20],[40,50],[80,90]])

#	-------------------------	Load the fits file	---------------------------

print("Reading "+fitsname)
necub,dkpc,theta0,phi0	=	fitld(fitsname,3.2)
print("Ne cube dimensions ")
print(necub.shape)
print("Spatial resolutions (kpc)")
print(dkpc)
print("Orientation (deg)")
print(np.rad2deg(theta0),np.rad2deg(phi0))

#	-------------------------	Execute tasks	-------------------------------

if (exmode=='profile'):
	print("\nGenerating radial electron density profiles...\n")
	cubene	= neprofinc(necub,dkpc,1.0,theta0,phi0,1.0,1.0,1.0)
	plot_nerad(cubene, incranges)

elif (exmode=='losdm'):
	print("\nEstimating LoS DMs...\n")
	losdms(fitsname,necub,dkpc,theta0,phi0,nfixpts,1.0,1.0,1.0)

elif (exmode=='pltdm'):
	print("\nPloting LoS DMs...\n")
	plotdms(fitsname,nfixpts,1.0,1.0,1.0,scalekpc)

elif (exmode=='dmscat'):
	print("\nPloting LoS DMs...\n")
	plotdm2d(fitsname,nfixpts,1.0,1.0,1.0,scalekpc)

else:
	print("\nHmm...What mode is that again...?\n")

plt.show()













































