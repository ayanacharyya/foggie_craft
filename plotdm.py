#
#	Plotting functions
#
#								AB, August 2024
#
#	Function list
#
#
#	--------------------------	Import modules	---------------------------

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import FuncFormatter
import matplotlib.colors as mpc
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit
import pickle as pkl
from collections import namedtuple
from globalpars import *
from nefns import *

mpl.rcParams['pdf.fonttype']	= 42
mpl.rcParams['ps.fonttype'] 	= 42
mpl.rcParams['savefig.dpi'] 	= 600
mpl.rcParams['font.family'] 	= 'sans-serif'
mpl.rcParams['font.size']		= 8

#	ne profiles for a given cube -- inclinatio, radius, density
#neradial		=	namedtuple('neradial',['logsm','logsfr','redshift','theta0','phi0','radkpc','theta','phi','inclination','neincrad'])

#	----------------------------------------------------------------------------------------------------------

def plot_nerad(radne,inclims):
	#	Plot Radial profile of ne within the inclination range
	
	fig 	= plt.figure(figsize=(6.0,4.5))
	
	colist	= ['b','r','k']
	shlist	= ['aqua','coral','grey']
	
	for i in range(0,len(inclims)):
		inclim	= inclims[i]
		relinds	= np.where((radne.inclination - np.deg2rad(inclim[0]))*(radne.inclination - np.deg2rad(inclim[1])) <= 0.0)	
		relincs	= radne.inclination[relinds]
		relne	= radne.neincrad[relinds]	
		medne	= np.nanmedian(relne, axis=0)
		meane	= np.nanmean(relne, axis=0)
		per16	= np.percentile(relne, 16.0, axis=0)
		per84	= np.percentile(relne, 84.0, axis=0)		
		plt.fill_between(radne.radkpc, per16, per84, color=shlist[i],alpha=0.2)
		plt.plot(radne.radkpc, medne, c=colist[i], ls='-', markersize=2, label=str(inclim[0])+"$^{\circ}$ < $i$ < "+str(inclim[1])+"$^{\circ}$")
		plt.plot(radne.radkpc, meane, c=colist[i], ls='--', markersize=2, label=str(inclim[0])+"$^{\circ}$ < $i$ < "+str(inclim[1])+"$^{\circ}$")
	
	plt.legend(loc="upper right",frameon=False)
	#plt.xlim([-1,205])
	plt.ylim([3.5e-7,1.5e-1])
	
	#plt.yscale('log')
	#plt.xscale('log')
	plt.xlabel('Radius (kpc)')
	plt.ylabel('$n_e$ (cm$^{-3}$)')
	plt.savefig("../plots/ne_radial_inc_"+str(inclim[0])+"_"+str(inclim[1])+".png",format='png',transparent=True)
		
	plt.show(block=False)
	
	return(0)
#	----------------------------------------------------------------------------------------------------------	



def plot_neimp(imfarr0, incarr0, losdm0, impbins, implim, fsize):
	#	Plot LoSDM vs impact factor for a given inclination range
	
	fig 	= plt.figure(figsize=(1.2*fsize,1.5*fsize))
	ax1	 	= fig.add_axes([0.15,0.15,0.84,0.30])
	ax2	 	= fig.add_axes([0.15,0.45,0.84,0.54])
	xtcarr	= [0.1,0.2,0.5,1,2,5,10,20,50,100]

	incbind	= []
	
	for ceninc in incvals:
		incindx	= incvals.index(ceninc)

		inincra	= np.where((incarr0-(ceninc - dinc/2.0))*(incarr0 - (ceninc+dinc/2.0)) < 0.0)
		imfarr	= imfarr0[inincra]
		losdm	= losdm0[inincra]

		print("Plotting DMs within inclination ",(ceninc - dinc/2.0), (ceninc + dinc/2.0))
	
		binned	= []
		impbineg= np.linspace(-1.0, np.log10(implim), impbins+1, dtype=float)
		dimp	= np.log10(implim)/impbins
		for i in range(0,impbins):
			bincen	= 10.0**((impbineg[i]+impbineg[i+1])/2)
			bindm	= losdm[(imfarr - 10.0**impbineg[i])*(imfarr - 10.0**impbineg[i+1]) < 0]
			medm	= np.nanmedian(bindm)
			binned.append([bincen, medm, np.nanmedian(np.abs(bindm-medm))])

		binned	= np.array(binned)
		#print(binned)
		
		#popt,pcov	= curve_fit(radialexp, binned[:,0], binned[:,1], sigma=binned[:,2], absolute_sigma=True)
		#perr 		= np.sqrt(np.diag(pcov))
		#print(popt, perr)

		ax2.fill_between(binned[:,0], binned[:,1]-binned[:,2], binned[:,1]+binned[:,2], color=shlist[incindx%len(shlist)],alpha=0.2)
		ax2.plot(binned[:,0], binned[:,1], clist[incindx%len(clist)]+lslist[incindx%len(lslist)], label=str((ceninc - dinc/2.0))+"$^{\circ}$ < $i$ < "+str((ceninc + dinc/2.0))+"$^{\circ}$")	
		#plt.plot(binned[:,0], radialexp(binned[:,0], *popt), clist[incindx%len(clist)]+':')

		ax1.plot(binned[:,0], binned[:,2]/binned[:,1], clist[incindx%len(clist)]+lslist[incindx%len(lslist)])
	
	plt.legend(loc="upper right",frameon=False)	
	ax2.set_yscale('log')
	ax1.set_xscale('log')
	ax2.set_xscale('log')
	ax2.set_ylim([0.8,1500])
	ax1.set_xticks(xtcarr)
	ax2.set_xticks(xtcarr)
	ax1.set_xticklabels(xtcarr)
	ax2.set_xticklabels([])
	#ax.set_yticks([1,10,100,1000])
	#ax.set_yticklabels([1,10,100,1000])
	ax1.set_xlabel('Impact parameter / $R_{eff}$')
	ax2.set_ylabel('Maximum DM (pc cm$^{-3}$)')
	ax1.set_ylabel('$\Delta$ DM / Maximum DM')
	#plt.savefig("../plots/ne_radial_inc_"+str(inclim[0])+"_"+str(inclim[1])+".png",format='png',transparent=True)
		
	plt.show(block=False)
	
	return(0)
#	----------------------------------------------------------------------------------------------------------	









































