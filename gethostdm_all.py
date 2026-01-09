#
#	Script for estimating FRB host galaxy DM from simulated electron density cubes
#
#								Originally by AB, August 2025
#                               Modified by AA, January 2026

#	--------------------------	Import modules	---------------------------

import os, sys, glob
import numpy as np
import pickle as pkl
from collections import namedtuple
from globalpars import *
from nefns import *
from plotdm import *
import subprocess
from mpi4py import MPI
from pathlib import Path
from datetime import timedelta, datetime
start_time = datetime.now()

# -------------------------------------------------------------------------------------------
def print_mpi(string):
    '''
    Function to print corresponding to each mpi thread
    '''
    comm = MPI.COMM_WORLD
    myprint_orig('[' + str(comm.rank) + '] {' + subprocess.check_output(['uname -n'],shell=True)[:-1].decode("utf-8") + '} ' + string + '\n')

# -------------------------------------------------------------------------------------------
def print_master(string):
    '''
    Function to print only if on the head node/thread
    '''
    comm = MPI.COMM_WORLD
    if comm.rank == 0: myprint_orig('[' + str(comm.rank) + '] ' + string + '\n')

# -------------------------------------------------------------------------------------------
def myprint_orig(text):
    '''
    Function to direct the print output to stdout or a file
    '''
    if not isinstance(text, list) and not text[-1] == '\n': text += '\n'
    if 'minutes' in text: text = fix_time_format(text, 'minutes')
    elif 'mins' in text: text = fix_time_format(text, 'mins')
    print(text)

# --------------------------------------------------------------------------------------------
def fix_time_format(text, keyword):
    '''
     Function to modify the way time is formatted in print statements
    '''
    arr = text.split(' ' + keyword)
    pre_time = ' '.join(arr[0].split(' ')[:-1])
    this_time = float(arr[0].split(' ')[-1])
    post_time = ' '.join(arr[1].split(' '))
    text = pre_time + ' %s' % (datetime.timedelta(minutes=this_time)) + post_time

    return text

# -------------------------------------------------------------------------------------------
def print_instructions():
    '''
	Print instructions to terminal
	'''
    print("\n            You probably need some assistance here!\n")
    print("\n Arguments are       --- <mode> <nfixpts> <extent/ckpc> <scale/kpc>\n")
    print(" Supported Modes are --- profile        (calculate electron density profiles)")
    print("                     --- losdm          (calculate LoS DMs)")
    print("                     --- pltdm          (Plot LoS DMs)")
    print("                     --- dmscat         (Plot LoS DMs)")
    print("\n            Now let's try again!\n")
	
    return(0)

# -----main code-----------------
if __name__ == '__main__':
    #	--------------------------	Read inputs	-------------------------------
    if(len(sys.argv)<3):
        print_instructions()
        sys.exit()

    exmode		=	sys.argv[1]					#	What to do	
    nfixpts     =   int(sys.argv[2])            #   Number of fixed points on each face to simulate LoSs
    extent     =   float(sys.argv[3])             #   Extent of datacube on either side of the center, in units of comoving kpc (to pick the right cubes from datadir)
    scalekpc	=	float(sys.argv[4])			#	Scale radius in kpc

    incranges	=	np.array([[0,20],[40,50],[80,90]])

    # --------domain decomposition; for mpi parallelisation-------------
    #list_of_fits = glob.glob(datadir + f'*El_number_density*{extent:.1f}ckpc*{scalekpc:.1f}kpc.fits') # all snapshots of this particular halo
    list_of_fits = glob.glob(losdir + f'*El_number_density*{extent:.1f}ckpc*{scalekpc:.1f}kpc*.npy') # all snapshots of this particular halo
    total_snaps = len(list_of_fits)

    comm = MPI.COMM_WORLD
    ncores = comm.size
    rank = comm.rank
    print_master(f'Total number of MPI ranks = {ncores} for total {total_snaps} snaps. ' + 'Starting at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.now()))
    comm.Barrier() # wait till all cores reached here and then resume

    split_at_cpu = total_snaps - ncores * int(total_snaps/ncores)
    nper_cpu1 = int(total_snaps / ncores)
    nper_cpu2 = nper_cpu1 + 1
    if rank < split_at_cpu:
        core_start = rank * nper_cpu2
        core_end = (rank+1) * nper_cpu2 - 1
    else:
        core_start = split_at_cpu * nper_cpu2 + (rank - split_at_cpu) * nper_cpu1
        core_end = split_at_cpu * nper_cpu2 + (rank - split_at_cpu + 1) * nper_cpu1 - 1

    # -------------loop over snapshots-----------------
    print_mpi('Operating on snapshots ' + str(core_start + 1) + ' to ' + str(core_end + 1) + ', i.e., ' + str(core_end - core_start + 1) + ' out of ' + str(total_snaps) + ' snapshots')
    start_index = 0

    for index in range(core_start + start_index, core_end + 1):
        start_time_this_snapshot = datetime.now()
        thisfile = Path(list_of_fits[index])
        fitsname = thisfile.stem
        if fitsname[-3:] == str(nfixpts): fitsname = fitsname[:-4]
        this_sim = fitsname.split('_')[:2]
        print_mpi('Doing snapshot ' + this_sim[0] + ' of halo ' + this_sim[1] + ' which is ' + str(index + 1 - core_start) + ' out of the total ' + str(core_end - core_start + 1) + ' snapshots...')

        #	-------------------------	Load the fits file	---------------------------
        if exmode in ['losdm', 'profile']:
            print_mpi("Reading "+fitsname)
            necub,dkpc,theta0,phi0	=	fitld(fitsname,3.2)
            print_mpi(f"Ne cube dimensions {necub.shape}")
            print_mpi(f"Spatial resolutions (kpc) {dkpc}")
            print_mpi(f"Orientation (deg) {np.rad2deg(theta0)},{np.rad2deg(phi0)}")

        #	-------------------------	Execute tasks	-------------------------------

        if (exmode=='profile'):
            print_mpi("\nGenerating radial electron density profiles...\n")
            cubene	= neprofinc(necub,dkpc,1.0,theta0,phi0,1.0,1.0,1.0)
            plot_nerad(cubene, incranges)

        elif (exmode=='losdm'):
            print_mpi("\nEstimating LoS DMs...\n")
            losdms(fitsname,necub,dkpc,theta0,phi0,nfixpts,1.0,1.0,1.0)

        elif (exmode=='pltdm'):
            print_mpi("\nPloting LoS DMs...\n")
            plotdms(fitsname,nfixpts,1.0,1.0,1.0,scalekpc)

        elif (exmode=='dmscat'):
            print_mpi("\nPloting LoS DMs...\n")
            plotdm2d(fitsname,nfixpts,1.0,1.0,1.0,scalekpc)

        else:
            print_mpi("\nHmm...What mode is that again...?\n")

        plt.show(block=False)
        print_mpi('This snapshots completed in %s' % timedelta(seconds=(datetime.now() - start_time_this_snapshot).seconds))

    # -----------------------------------------------------------------------------------
    if ncores > 1: print_master('Parallely: time taken for ' + str(total_snaps) + ' snapshots with ' + str(ncores) + ' cores was %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
    else: print_master('Serially: time taken for ' + str(total_snaps) + ' snapshots with ' + str(ncores) + ' core was %s' % timedelta(seconds=(datetime.now() - start_time).seconds))