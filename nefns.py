#
#	Functions for analysing electron density in simulated cubes
#
#								AB, August 2024
#
#	Function list
#
#	def fitld(fitsname):
#		Reads a FITS cube and returns a 3D numpy array and the spatial resolutions in kpc			
#
#	def neprofile(necube,dkpc,theta,phi):
#		Calculate electron density profiles along theta, phi
#	
#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pickle as pkl
from collections import namedtuple
from itertools import combinations
from globalpars import *
from auxfns import *
from plotdm import *

mpl.rcParams['pdf.fonttype']	= 42
mpl.rcParams['ps.fonttype'] 	= 42
mpl.rcParams['savefig.dpi'] 	= 600
mpl.rcParams['font.family'] 	= 'sans-serif'
mpl.rcParams['font.size']		= 8

#	----------------------------------------------------------------------------------------------------------
def fitld(fitsname, fsize):	
#	Reads a FITS cube and returns a 3D numpy array and the spatial resolutions in kpc
	
	fitsfile	=	fits.open(datadir+fitsname+".fits")
	fitshdr		=	fitsfile[0].header
	dxkpc		=	fitshdr['CDELT1']
	dykpc		=	fitshdr['CDELT2']
	dzkpc		=	fitshdr['CDELT3']
	nvecx		=	fitshdr['HIERARCH NORMAL_UNIT_VECTOR1']
	nvecy		=	fitshdr['HIERARCH NORMAL_UNIT_VECTOR2']
	nvecz		=	fitshdr['HIERARCH NORMAL_UNIT_VECTOR3']
	necube		=	fitsfile[0].data
	fitsfile.close()

	reskpc		=	(dxkpc,dykpc,dzkpc)
	theta0		=	np.arctan2(nvecz, np.sqrt(nvecx*nvecx + nvecy*nvecy))
	phi0		=	np.arctan2(nvecy, nvecx)
	
	fig 	= plt.figure(figsize=(3*fsize, fsize))
	ax0 		= fig.add_axes([0.05,0.10,0.31,0.88])
	ax0.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	plt.imshow(np.nansum(necube,axis=0),interpolation='none',origin='lower')
	
	ax1 		= fig.add_axes([0.36,0.10,0.31,0.88])
	ax1.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	plt.imshow(np.nansum(necube,axis=1),interpolation='none',origin='lower')
	ax1.set_yticks([])
	
	ax2 		= fig.add_axes([0.67,0.10,0.31,0.88])
	ax2.tick_params(axis="both",direction="in",bottom=True,right=True,top=True,left=True)
	plt.imshow(np.nansum(necube,axis=2),interpolation='none',origin='lower')
	ax2.set_yticks([])

	plt.show(block=False)
	
	return (necube,reskpc,theta0,phi0)
#	----------------------------------------------------------------------------------------------------------



def neprofinc(necube,dkpc,dangdeg,theta0,phi0,logsm,logsfr,redshift):	
#	Calculate radial electron density profile and inclinations along many radial directions 
	
	dangrad	= np.deg2rad(dangdeg)
	cenpx	= np.array(necube.shape)/2
	radius	= np.arange(0,np.amin(cenpx),np.amax(dkpc),dtype=np.float32)
	thetas	= np.linspace(-(np.pi/2.0), (np.pi/2.0), int(np.rint(np.pi/dangrad)), dtype=np.float32)
	angls	= []
	for theta in thetas:
		phis	= np.arange(0.0, 2*np.pi*np.abs(np.cos(theta))+dangrad/2.0, dangrad, dtype=np.float32)
		incs	= np.arccos(np.sin(theta)*np.sin(theta0) + np.cos(theta)*np.cos(theta0)*np.cos(phis-phi0))
		incs	= np.arcsin(np.sin(incs))
		for i in range(0,len(phis)):
			angls.append([theta,phis[i],incs[i]])			
	angls	= np.array(angls)
	
	nepros	= np.zeros((len(angls), len(radius)), dtype=np.float32)
	
	for i in range(0,len(angls)):
		nepros[i]	= neprofile(necube,dkpc,angls[i,0],angls[i,1],cenpx,radius)	
	
	neres	= neradial(logsm,logsfr,redshift,theta0,phi0,radius,angls[:,0],angls[:,1],angls[:,2],nepros)
		
	return(neres)
#	------------------------------------------------------------------------------------------------------



def intnelos(necube,dkpc,xyz1,xyz2):	
#	Integrate ne along a given LoS
#	Arguments:	ne cube
#				cell/pixel size in kpc
#	
	vec		= xyz2 - xyz1
	lenlos	= np.sqrt(np.sum(vec**2))		
	#	Parameter of the parametric equation of a straight line
	tarr 	= np.arange(0.0, 1.0, 1.0/lenlos)
		
	xarr	= np.rint(vec[0]*tarr + xyz1[0] + necube.shape[0]/2).astype(int) % necube.shape[0]
	yarr	= np.rint(vec[1]*tarr + xyz1[1] + necube.shape[1]/2).astype(int) % necube.shape[1]
	zarr	= np.rint(vec[2]*tarr + xyz1[2] + necube.shape[2]/2).astype(int) % necube.shape[2]
	
	#print(np.array([xarr,yarr,zarr]).T)
	dmlos	= np.nansum(necube[xarr, yarr, zarr])*np.mean(dkpc)*1.0e3

	return (dmlos)
#	------------------------------------------------------------------------------------------------------



def losdms(fitsname,necube,dkpc,theta0,phi0,nfixpts,logsm,logsfr,redshift):	
#	Calculate DMs integrating over LoS

	dmarr	= []

	ptpairs	= np.random.uniform(-necube.shape[0]/2, necube.shape[0]/2, size=(6,2,nfixpts))
	
	ptsxl	= np.array([-(necube.shape[0]/2)*np.ones(nfixpts,dtype=int), ptpairs[0,0], ptpairs[0,1]]).T
	ptsxr	= np.array([ (necube.shape[0]/2)*np.ones(nfixpts,dtype=int), ptpairs[1,0], ptpairs[1,1]]).T
	ptsyl	= np.array([ptpairs[2,0], -(necube.shape[0]/2)*np.ones(nfixpts,dtype=int), ptpairs[2,1]]).T
	ptsyr	= np.array([ptpairs[3,0],  (necube.shape[0]/2)*np.ones(nfixpts,dtype=int), ptpairs[3,1]]).T
	ptszl	= np.array([ptpairs[4,0], ptpairs[4,1], -(necube.shape[0]/2)*np.ones(nfixpts,dtype=int)]).T
	ptszr	= np.array([ptpairs[5,0], ptpairs[5,1],  (necube.shape[0]/2)*np.ones(nfixpts,dtype=int)]).T
	
	ptslist	= [ptsxl, ptsxr, ptsyl, ptsyr, ptszl, ptszr]
	plcom	= combinations(ptslist, 2)
	
	for ac in plcom:
		for pt1 in ac[0]:
			for pt2 in ac[1]:
				vec		= pt2 - pt1
				incdeg	= inclinvec(vec, theta0, phi0)
				impf 	= impactfac(pt1, pt2)*np.mean(dkpc)
				dmlos	= intnelos(necube,dkpc,pt1,pt2)
				#print(vec, incdeg, impf, dmlos)
				dmarr.append([incdeg, impf, dmlos])
    
	dmarr	= np.array(dmarr)
	
	print("Total number of LoS = ",dmarr.shape[0])
	print("Saving LoS DMs to ",fitsname)
	np.save(losdir+fitsname+"_"+str(nfixpts)+".npy",dmarr)
	
	return(0)
#	------------------------------------------------------------------------------------------------------



def plotdms(fitsname,nfixpts,logsm,logsfr,redshift,scalekpc):	
#	Plots maximum DMs along different LoSs
	
	dmarr	= np.load(losdir+fitsname+"_"+str(nfixpts)+".npy")	
	print("Total number of LoS = ",dmarr.shape[0])	
	
	plot_neimp(dmarr[:,1]/scalekpc, dmarr[:,0], dmarr[:,2], impbins, maximpa, 3.2)

	return(0)
#	------------------------------------------------------------------------------------------------------








































































