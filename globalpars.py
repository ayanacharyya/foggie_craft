#	Globally defined parameters which should not be frequently modified
import numpy as np
from collections import namedtuple

#root_dir    =   "../data_260310/"
#root_dir = "../"
root_dir = "/Users/acharyya/Library/CloudStorage/GoogleDrive-ayan.acharyya@inaf.it/My Drive/FOGGIE_CRAFT/"
#root_dir = "/nobackupp19/aachary2/foggie_craft/"

#datadir     =   root_dir                    #   Location of the FITS cubes
datadir     =   root_dir + "data/"                #   Location of the FITS cubes
losdir      =   root_dir + "losdms/"        #   Location of the LoS DMs
radialdir   =   datadir + "radial_profiles/"                 #   Location of the radial profiles
plotdir     =   root_dir + "plots/"        #   Location of the LoS DMs
plotradial  =   plotdir + "density_profiles/"                   # 

impbinegs   =   np.array([0,2,4,8,16,32,64,128,256])              #   Impact parameter bins
#impbinegs   =   np.array([0,1,2,4,8,16,32,64,128,256])                #   Impact parameter bins in units of r_eff
maxdmcol    =   205.0                                             #   Maximum DM for colour scale
impbins     =   20                  #   Number of bins in impact factor
maximpa     =   100.0               #   Maximum value of impact factor in units of R_eff

incvals     =   [5.0,45.0,85.0]                                   #   Central values of inclination bins in deg
dinc        =   10.0                                              #   Width of the inclination bins in deg
radbins     =   [0, 1, 2, 4, 8, 16, 32, 64, 128]                #   Radial )bin edges (in kpc)
dm_ticks    =   [1, 3, 10, 30, 100, 300]

clist       =   ['b', 'r', 'k']
shlist      =   ['aqua', 'lightsalmon', 'lightgrey']
lslist      =   ['-', '--', ':']
marklist    =   ['*', 's', 'o']

scale_fit_thresh = 3.0                                      # threshold for scaling relation robust fitting, in sigma

# ----------------global variables for plotting routines-----------------
#	A named tuple to store various parameters related to radial ne profile
neradial		=	namedtuple('neradial',['logsm','logsfr','redshift','theta0','phi0','radkpc','theta','phi','inclination','neincrad'])

all_lsm_bin_edges = [8.5,8.75,9.0,9.25,9.5,9.75,10.0,10.25,10.5,10.75,11.0,11.25,11.5]
all_lsfr_bin_edges	=[-3.0,0.0,0.5,1.0,1.5,2.0]





















