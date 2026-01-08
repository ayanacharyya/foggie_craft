import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pickle as pkl
from collections import namedtuple

#	Globally defined parameters which should not be frequently modified

#root_dir = "../"
#root_dir = "/Users/acharyya/Library/CloudStorage/GoogleDrive-ayan.acharyya@inaf.it/My Drive/FOGGIE_CRAFT/"
root_dir = "/nobackupp19/aachary2/foggie_craft/"
datadir     =   root_dir + "data/"     #   Location of the FITS cubes
losdir      =   root_dir + "losdms/"        #   Location of the LoS DMs
impbins     =   20                  #   Number of bins in impact factor
maximpa     =   100.0               #   Maximum value of impact factor in units of R_eff
incvals     =   [15.0,45.0,75.0]    #   Central values of inclination bins in deg
dinc        =   30.0                #   Width of the inclination bins in deg
clist       =   ['b', 'r', 'k']
shlist      =   ['aqua', 'lightsalmon', 'lightgrey']
lslist      =   ['-', '-', '-']
























