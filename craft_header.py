"""

    Title :      header
    Notes :      Header file for importing packages/modules and declaring global variables required for working with FOGGIE-CRAFT code.
    Author :     Ayan Acharyya
    Started :    31-03-26

"""

import os, sys, argparse, re, subprocess, time, math, glob
import numpy as np
import multiprocessing as mproc
from datetime import timedelta, datetime

import pickle as pkl
from collections import namedtuple

from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.path import Path as mpl_Path
from matplotlib import cm as mpl_cm
from matplotlib import colors as mplcolors
from matplotlib import ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import mplcyberpunk

from pathlib import Path
from importlib import reload

#from mpi4py import MPI

from scipy import optimize as op
from scipy import stats
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy.optimize import curve_fit, fminbound
from scipy.ndimage import gaussian_filter
from scipy.stats import binned_statistic
from scipy.stats import binned_statistic_2d

from astropy.io import ascii, fits
from astropy.table import Table
from astropy.stats import gaussian_fwhm_to_sigma as gf2s
from astropy import convolution as con
from astropy import units as u

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)

from foggie_craft_utils.get_run_loc_etc import *

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

