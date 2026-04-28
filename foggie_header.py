#!/usr/bin/env python3

"""

    Title :      header
    Notes :      Header file for importing packages/modules and declaring global variables required for working with FOGGIE code.
    Author :     Ayan Acharyya
    Started :    January 2021

"""

import numpy as np
import multiprocessing as mproc
import seaborn as sns
import os
import sys
import argparse
import re
import subprocess
import time
import math
import shutil
import copy
import glob
import random
import collections, itertools

import pickle as pkl
from collections import namedtuple

import subprocess

from matplotlib import pyplot as plt
#plt.style.use('seaborn-whitegrid')
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['xtick.top'] = True

from matplotlib import colors as mplcolors
from matplotlib import patheffects as fx
from matplotlib.colors import LogNorm
from matplotlib import image as mpimg
from matplotlib.path import Path as mpl_Path
from matplotlib import cm as mpl_cm
import mplcyberpunk


from pathlib import Path
from importlib import reload

from mpi4py import MPI

from numpy import exp
from scipy import optimize as op
from scipy import stats
from scipy.interpolate import interp1d
from scipy.interpolate import RegularGridInterpolator as RGI
from scipy.interpolate import LinearNDInterpolator as LND
from scipy.special import erf
from scipy.optimize import curve_fit, fminbound
from scipy.ndimage import gaussian_filter
from scipy.stats import binned_statistic

from astropy.io import ascii, fits
from astropy.table import Table
from astropy.stats import gaussian_fwhm_to_sigma as gf2s
from astropy import convolution as con
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM, Planck13, z_at_value

from operator import itemgetter
from collections import defaultdict
import cmasher as cmr
from uncertainties import ufloat, unumpy

import datashader as dsh
from datashader.utils import export_image
from datashader import transfer_functions as dstf

datashader_ver = float(dsh.__version__.split('.')[1])
if datashader_ver > 11: from datashader.mpl_ext import dsshow

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)

import yt
from yt.units import *

from foggie_craft_utils.get_run_loc_etc import *
from foggie_craft_utils.consistency import *
from foggie_craft_utils.yt_fields import *
from foggie_craft_utils.foggie_load import *
from foggie_craft_utils.get_proper_box_size import get_proper_box_size
from foggie_craft_utils.util import *

from datetime import timedelta, datetime

# ------------declaring constants to be used globally-----------
c = 3e5  # km/s
H0 = 70.  # km/s/Mpc Hubble's constant
planck = 6.626e-27  # ergs.sec Planck's constant
nu = 5e14  # Hz H-alpha frequency to compute photon energy approximately
Mpc_to_m = 3.08e22
Mpc_to_cm = Mpc_to_m * 100
kpc_to_cm = Mpc_to_cm / 1000

#alpha_B = 3.46e-19  # m^3/s OR 3.46e-13 cc/s, Krumholz & Matzner (2009) for 7e3 K
alpha_B = 2.59e-19  # m^3/s OR 2.59e-13 cc/s, for Te = 1e4 K, referee quoted this values
k_B = 1.38e-23  # m^2kg/s^2/K
G = 6.67e-11  # Nm^2/kg^2
eps = 2.176e-18  # Joules or 13.6 eV
m_H = 1.67e-27  # kg; mass of proton

# ------------declaring overall paths (can be modified on a machine/user basis)-----------
HOME = os.getenv('HOME')
try:
    if not os.path.exists(HOME+'/Work/astro/ayan_codes'): # if the code directory does not exist in current home, then it must exist in /pleiades home
        HOME = '/pleiades/u/' + os.getenv('USER')
except:
    pass

mappings_lab_dir = HOME + '/Mappings/lab/'  # if you are producing the MAPPINGS grid,
# this is where your MAPPINGS executable .map51 is installed,
# otherwise, this is where your MAPPINGS grid and your emission line list is
mappings_input_dir = HOME + '/Mappings/HIIGrid306/Q/inputs/'  # if you are producing the MAPPINGS grid,
# this is where your MAPPINGS input/ directory is
# otherwise, ignore this variable
sb99_dir = HOME + '/SB99-v8-02/output/'  # this is where your Starburst99 model outputs reside
# this path is used only when you are using compute_hiir_radii.py or lookup_flux.py
sb99_model = 'starburst11'  # for fixed stellar mass input spectra = 1e6 Msun, run up to 10 Myr
sb99_mass = 1e6  # Msun, mass of star cluster in given SB99 model

# ------------declaring list of ALL simulations present locally (on HD)-----------
all_sims_dict = {'8508': [('8508', 'RD0042'), ('8508', 'RD0039'), ('8508', 'RD0031'), ('8508', 'RD0030'), ('8508', 'DD2288'), ('8508', 'DD2289')], \
                 '5036': [('5036', 'RD0039'), ('5036', 'RD0031'), ('5036', 'RD0030'), ('5036', 'RD0020')], \
                 '5016': [('5016', 'RD0042'), ('5016', 'RD0039'), ('5016', 'RD0031'), ('5016', 'RD0030'), ('5016', 'RD0020')], \
                 '4123': [('4123', 'RD0031'), ('4123', 'RD0030')], \
                 '2878': [('2878', 'RD0020'), ('2878', 'RD0018')], \
                 '2392': [('2392', 'RD0030')], \
    } # all snapshots in the HD

projection_dict = {'x': ('y', 'z', 'x'), 'y':('z', 'x', 'y'), 'z':('x', 'y', 'z')} # which axes are projected for which line of sight args.projection

# -----------declaring/modifying colormaps to be ued for certain properties throughout my code------------
# individually comment out following lines to keep the original color_map as defined in foggie.utils.consistency

#density_color_map = 'viridis'

velocity_discrete_cmap = 'coolwarm'

temperature_color_list = ("darkred", "#d73027", "darkorange", "#ffe34d")
temperature_color_map = sns.blend_palette(temperature_color_list, as_cmap=True)

metal_color_list = ("#4575b4", "#984ea3", "#984ea3", "#d73027", "darkorange", "#ffe34d")
#metal_color_list = ("black", "#4575b4", "#984ea3", "#d73027", "darkorange", "#ffe34d")
metal_color_map = sns.blend_palette(metal_color_list, as_cmap=True)
#metal_color_map = 'viridis'

metal_colors_mw = sns.blend_palette(metal_color_list, n_colors=6)
metal_discrete_cmap_mw = mplcolors.ListedColormap(metal_colors_mw)
metal_color_key_mw = collections.OrderedDict()
for i in np.arange(np.size(metal_color_labels_mw)): metal_color_key_mw[metal_color_labels_mw[i]] = to_hex(metal_colors_mw[i])

# --- Enable autoreload when running in IPython ---
try:
    from IPython import get_ipython
    ipython = get_ipython()
    if ipython is not None:
        ipython.run_line_magic('load_ext', 'autoreload')
        ipython.run_line_magic('reload_ext', 'autoreload')
        ipython.run_line_magic('autoreload', '2')
except Exception:
    pass

#-------- set variables and dictionaries such that they are available to other scripts importing this script-----------
field_dict = {'radius':('gas', 'radius_corrected'), 'density':('gas', 'density'), 'mass':('gas', 'mass'), \
              'metal':('gas', 'metallicity'), 'temp':('gas', 'temperature'), 'vrad':('gas', 'radial_velocity_corrected'), 'vdisp': ('gas', 'velocity_dispersion_3d'), \
              'phi_L':('gas', 'angular_momentum_phi'), 'theta_L':('gas', 'angular_momentum_theta'), 'volume':('gas', 'volume'), \
              'phi_disk': ('gas', 'phi_pos_disk'), 'theta_disk': ('gas', 'theta_pos_disk'), 'el_density':('gas', 'El_number_density')}
if yt.__version__[0]=='3':
    field_dict['mass'] = ('gas','cell_mass')
    field_dict['volume'] = ('gas', 'cell_volume')
unit_dict = {'radius':'kpc', 'rad_re':'', 'density':'g/cm**3', 'metal':r'Zsun', 'temp':'K', 'vrad':'km/s', 'phi_L':'deg', 'theta_L':'deg', 'PDF':'', 'mass':'Msun', 'stars_mass':'Msun', 'ystars_mass':'Msun', 'ystars_age':'Gyr', 'gas_frac':'', 'gas_time':'Gyr', 'volume':'pc**3', 'phi_disk':'deg', 'theta_disk':'deg', 'vdisp':'km/s', 'el_density':'cm**-3'}
labels_dict = {'radius':'Radius', 'rad_re':'Radius/R_e', 'density':'Density', 'metal':'Metallicity', 'temp':'Temperature', 'vrad':'Radial velocity', 'phi_L':r'$\phi_L$', 'theta_L':r'$\theta_L$', 'PDF':'PDF', 'gas_frac':'Gas fraction', 'gas_time':'Gas consumption timescale', 'phi_disk':'Azimuthal Angle', 'theta_disk':r'$\theta_{\mathrm{diskrel}}$', 'vdisp':'Gas velocity dispersion', 'el_density': 'Electron number density'}
